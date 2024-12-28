//  monalisa
//
//  Created by Paul O. Lewis 12/20/2024.
//  Inspired by app created by Roger Johansson:
//     https://rogerjohansson.blog/2008/12/07/
//     genetic-programming-evolution-of-mona-lisa/
//
//  Change imagefile in monalisa.conf file to specify reference image
//  400x400 image is good; larger images will take longer per generation
//  Press Q key on any window to quit after current generation
//
//  Simulates natural selection in a "population" the individuals of which
//  have a genome that is an image of the same size as the reference image.
//  Each generation, the best nreprod of the nindivs individuals in the
//  population are preserved and the remaining nindivs - nreprod individuals
//  are replaced by copying from randomly-chosen preserved individuals.
//  All individuals undergo mutation, which changes the color of a single
//  pixel in their genome.
//
//  While this program allows multithreading, you will probably find it
//  does not help because the "fitness" calculation is so cheap that there
//  is more overhead than processing for the multithreaded version.
//
//  This program depends on the following libraries:
//    boost (program_options):
//      https://www.boost.org/doc/libs/1_87_0/doc/html/program_options.html
//    CImg: https://cimg.eu/index.html
//  It also requires X11 for the graphical display (e.g. XQuartz for MacOS)
//  A makefile is provided that can be adjusted for your system.

#include <fstream>
#include <limits>
#include <vector>
#include <string>
#include <cmath>
#include <thread>
using namespace std;

string program_name        = "monalisa";
unsigned major_version     = 1;
unsigned minor_version     = 0;
string author              = "Paul O. Lewis";
string minor_version_date  = "20-Dec-2024";

#include <boost/format.hpp>
#include <boost/program_options.hpp>
using boost::format;
using boost::program_options::variables_map;
using boost::program_options::options_description;
using boost::program_options::value;
using boost::program_options::store;
using boost::program_options::parse_command_line;
using boost::program_options::parsed_options;
using boost::program_options::parse_config_file;
using boost::program_options::reading_file;
using boost::program_options::notify;

#include <boost/random/gamma_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include "xproj.hpp"
#include "stopwatch.hpp"
#include "lot.hpp"
using namespace proj;

Lot lot;

#define cimg_display 1
#define cimg_use_png 0
#include "CImg.h"
using namespace cimg_library;

typedef unsigned char img_t;

CImg<img_t> saved_image = CImg<img_t>::empty();

struct ThreadSchedule {
    ThreadSchedule(unsigned nt, unsigned ni) : _nthreads(nt), _nindivs(ni) {
        // first is index of first individual in range
        // last is index of last individual in range plus 1
        assert(_nthreads > 0);
        assert(_nthreads <= _nindivs);
        _first.resize(_nthreads);
        _last.resize(_nthreads);
        if (_nthreads == 1) {
            _first[0] = 0;
            _last[0] = _nindivs;
        }
        else {
            // num_threads > 1
            unsigned stride = _nindivs/_nthreads;
            unsigned r = _nindivs % _nthreads;
            unsigned first = 0;
            unsigned last = first;
            for (unsigned i = 0; i < _nthreads; i++) {
                first = last;
                last = first + stride;
                if (r > 0) {
                    last++;
                    r--;
                }
                if (last > _nindivs)
                    last = _nindivs;
                _first[i] = first;
                _last[i] = last;
            }
        }
    }
    
    string summary() {
        string s = "";
        if (_nthreads == 1) {
            s += "All individuals handled by the single thread\n";
        }
        else {
            s += "First and last individual handled by each thread\n";
            s += str(format("  %12s %12s %12s %12s\n") % "thread" % "first" % "last" % "n");
            unsigned nsum = 0;
            for (unsigned i = 0; i < _nthreads; i++) {
                nsum += _last[i] - _first[i];
                s += str(format("  %12d %12d %12d %12d\n") % (i+1) % (_first[i] + 1) % _last[i] % (_last[i] - _first[i]));
            }
            s += str(format("  %12s %12s %12s %12d\n") % " " % " " % " " % nsum);
        }
        return s;
    }
    
    unsigned _nthreads;
    unsigned _nindivs;
    vector<unsigned> _first;
    vector<unsigned> _last;
};

struct Operators {
    Operators() : _sum_weights(0.0) {
    }
    
    void add(unsigned index, float weight) {
        assert(_weights.count(index) == 0);
        _weights[index] = weight;
    }
    
    float getSumWeights() {
        return _sum_weights;
    }
    
    void finalize() {
        _sum_weights = 0.0;
        for (auto w : _weights) {
            _sum_weights += w.second;
        }
    }
    
    unsigned choose(float u) {
        assert(_sum_weights > 0.0);
        float cumprob = 0.0;
        unsigned index = 0;
        for (auto w : _weights) {
            float prob = w.second/_sum_weights;
            cumprob += prob;
            if (u < cumprob) {
                index = w.first;
                break;
            }
        }
        assert(index > 0);
        return index;
    }
    
    map<unsigned,float> _weights;
    float _sum_weights;
};

class Individual {
    public:
        Individual() {
            clear();
        }
        
        ~Individual() {
            clear();
        }
        
        void clear() {
            _score          = numeric_limits<float>::max();
            _prev_score     = _score;
            _genome         = CImg<img_t>::empty();
            _x              = 0;
            _y              = 0;
            _prev_red       = 255;
            _prev_green     = 255;
            _prev_blue      = 255;
            _curr_red       = 255;
            _curr_green     = 255;
            _curr_blue      = 255;
        }
        
        void initialize(const CImg<img_t> & ref) {
            // Make a copy of ref
            _genome = ref;
            
            // Initialize colors
            img_t * r  = _genome.data(0, 0, 0, 0);
            img_t * g  = _genome.data(0, 0, 0, 1);
            img_t * b  = _genome.data(0, 0, 0, 2);
            for(int row = 0; row < _genome.height(); row++) {
                for(int col = 0; col < _genome.width(); col++) {
                    *r = (img_t)max_rgb*::lot.uniform();
                    *g = (img_t)max_rgb*::lot.uniform();
                    *b = (img_t)max_rgb*::lot.uniform();

                    ++r;
                    ++g;
                    ++b;
                }
            }

            r  = _genome.data(_x, _y, 0, 0);
            g  = _genome.data(_x, _y, 0, 1);
            b  = _genome.data(_x, _y, 0, 2);
            _prev_red   = *r;
            _prev_green = *g;
            _prev_blue  = *b;
            _curr_red   = *r;
            _curr_green = *g;
            _curr_blue  = *b;
        }
        
        void mutateColor(img_t & color) {
            float u = (float)::lot.uniform();
            float offset = (float)(-0.5*delta_rgb + delta_rgb*u);
            float value = color + offset;
            if (value < 0.0 || value > max_rgb) {
                if (value < 0.0)
                    value = -value;
                if (value > max_rgb)
                    value = max_rgb - (value - max_rgb);
            }
            //TODO: value can never exactly equal max_rgb
            color = (img_t)value;
        }
        
        float getScore() {return _score;}
        
        void mutate() {
            assert(_curr_red == _prev_red);
            assert(_curr_green == _prev_green);
            assert(_curr_blue == _prev_blue);
            
            // Choose a row and column at random to mutate
            _x = ::lot.randint(0, refh - 1);
            _y = ::lot.randint(0, refw - 1);
                        
            // Mutate all three channels for pixel at _x,_y
            img_t * r  = _genome.data(_x, _y, 0, 0);
            img_t * g  = _genome.data(_x, _y, 0, 1);
            img_t * b  = _genome.data(_x, _y, 0, 2);
            _prev_red   = *r;
            _prev_green = *g;
            _prev_blue  = *b;
            mutateColor(*r);
            mutateColor(*g);
            mutateColor(*b);
            _curr_red   = *r;
            _curr_green = *g;
            _curr_blue  = *b;
        }
        
        void replaceGenome(CImg<img_t> & replacement) {
            _genome = replacement;
        }
                
        void display(CImgDisplay & disp, bool start_message) {
            if (start_message) {
                // Save _genome so that it can be replaced when the message is supposed to go away
                ::saved_image = _genome;
                
                // Draw the text on an empty image in order to determine the rendered dimensions
                CImg<unsigned char> dummy;
                vector<img_t> fg = {255, 255, 255};
                vector<img_t> bg = {0, 0, 0};
                dummy.draw_text(0, 0, "Press S to begin", &fg[0], &bg[0], 1, 23);
                unsigned w = dummy.width();
                unsigned h = dummy.height();
                
                unsigned x = _genome.width()/2.0 - w/2.0;
                unsigned y = _genome.height()/2.0 - h/2.0;
                _genome.draw_text(x, y, "Press S to begin", &fg[0], &bg[0], 1.0, 23);
            }
            disp = _genome;
        }
        
        void saveToFile(string fn, int version_number) {
            _genome.save(fn.c_str(), version_number);
        }
        
        void resetPrevScore() {
            _prev_score = _score;
            _prev_red   = _curr_red;
            _prev_green = _curr_green;
            _prev_blue  = _curr_blue;
        }
        
        double calcScoreXY(const CImg<img_t> & refimg) {
            // Assume all pixels are the same except pixel at (_x,_y)
            float score = pow(_prev_score,2.0);
            const img_t * r0 = refimg.data(_x, _y, 0, 0);
            const img_t * g0 = refimg.data(_x, _y, 0, 1);
            const img_t * b0 = refimg.data(_x, _y, 0, 2);
            float srprev = pow(_prev_red   - *r0, 2.0);
            float sgprev = pow(_prev_green - *g0, 2.0);
            float sbprev = pow(_prev_blue  - *b0, 2.0);
            score -= (srprev + sgprev + sbprev);
            const img_t * r  = _genome.data(_x, _y, 0, 0);
            const img_t * g  = _genome.data(_x, _y, 0, 1);
            const img_t * b  = _genome.data(_x, _y, 0, 2);
            float sr = pow(*r - *r0, 2.0);
            float sg = pow(*g - *g0, 2.0);
            float sb = pow(*b - *b0, 2.0);
            score += (sr + sg + sb);
            _score = sqrt(score);

            return _score;
        }

        double calcScore(const CImg<img_t> & refimg) {
            // Fully compute score
            float score = 0.0f;
            const img_t * r0 = refimg.data(0, 0, 0, 0);
            const img_t * g0 = refimg.data(0, 0, 0, 1);
            const img_t * b0 = refimg.data(0, 0, 0, 2);
            const img_t * r  = _genome.data(0, 0, 0, 0);
            const img_t * g  = _genome.data(0, 0, 0, 1);
            const img_t * b  = _genome.data(0, 0, 0, 2);
            for(int row = 0; row < _genome.height(); row++) {
                for(int col = 0; col < _genome.width(); col++) {
                    float sr = pow(*r - *r0, 2.0);
                    float sg = pow(*g - *g0, 2.0);
                    float sb = pow(*b - *b0, 2.0);
                    score += sr + sg + sb;
                    ++r0; ++r;
                    ++g0; ++g;
                    ++b0; ++b;
                }
            }
            _score = sqrt(score);
            return _score;
        }

        static void calcCumProb() {
            operators.add(OP_MUTATE, mutate_weight);
            operators.finalize();
        }

        static unsigned refw;
        static unsigned refh;

        static float max_rgb;
        static float delta_rgb;

        static float mutate_weight;

        static bool  untouchable;
        
    private:
    
        enum OPERATOR {
            OP_MUTATE    = 1,
        };
        
        static Operators operators;

        CImg<img_t> _genome;
        
        float       _score;
        float       _prev_score;
        
        unsigned    _x;
        unsigned    _y;
        
        img_t       _prev_red;
        img_t       _prev_green;
        img_t       _prev_blue;
        
        img_t       _curr_red;
        img_t       _curr_green;
        img_t       _curr_blue;
};

void initializePopulation(vector<Individual> & pop, const CImg<img_t> & ref) {
    for (auto it = pop.begin(); it != pop.end(); ++it) {
        it->initialize(ref);
    }
}

// Default settings
string reference_image_filename                 = "";
string output_filename                          = "output.txt";
bool verbose                                    = false;
bool display_all_time_best                      = true;
unsigned report_every                           = 100;
unsigned save_every                             = 1000;
unsigned stop_at_generation                     = numeric_limits<unsigned>::max();
unsigned nthreads                               = 1;
unsigned rnseed                                 = 0;
unsigned nindivs                                = 20;
unsigned nreprod                                = 10;

float    Individual::delta_rgb                  = 100;

float    Individual::mutate_weight              = 0.0;

float    Individual::max_rgb                    = 255;
bool     Individual::untouchable                = true;
unsigned Individual::refw                       = 0;
unsigned Individual::refh                       = 0;

Operators Individual::operators;

vector< pair<unsigned, float> > progress;

void processCommandLineOptions(int argc, const char * argv[]) {
    variables_map vm;
    options_description desc("Allowed options");
    desc.add_options()
    ("help,h", "produce help message")
    ("version,v", "show program version")
    ("imagefile",  value(&reference_image_filename), "name of a JPG image to use as a reference")
    ("outfile",  value(&output_filename), "filename in which to save console output")
    ("verbose",  value(&verbose)->default_value(false), "verbose yields more output")
    ("alltimebest", value(&display_all_time_best)->default_value(true), "if true, best individual ever seen is displayed; if false, best individual from last generation is displayed")
    ("reportevery",  value(&report_every)->default_value(100), "how many generations to skip before reporting progress")
    ("saveevery",  value(&save_every)->default_value(1000), "how many generations to skip before saving best individual")
    ("stopat",  value(&stop_at_generation)->default_value(numeric_limits<unsigned>::max()), "stop after this many generations")
    ("nthreads",  value(&nthreads)->default_value(1), "number of threads to use in parallelizing fitness calculation ")
    ("rnseed",  value(&rnseed)->default_value(0), "pseudorandom number generator seed")
    ("nindivs",  value(&nindivs)->default_value(20), "number of individuals in the population")
    ("nreprod",  value(&nreprod)->default_value(10), "number of individuals allowed to reproduce (must be less than or equal to nindivs)")
    ("deltargb",  value(&Individual::delta_rgb)->default_value(100), "window size for proposals to change one color channel")
    ("untouchable",  value(&Individual::untouchable)->default_value(true), "best individual is preserved each generation and not allowed to be touched by mutation")
    ("wtmutate",  value(&Individual::mutate_weight)->default_value(0.0), "weight for proposals to change the color of a pixel")
    ;
    
    store(parse_command_line(argc, argv, desc), vm);
    try {
        const parsed_options & parsed = parse_config_file< char >("monalisa.conf", desc, false);
        store(parsed, vm);
    }
    catch(reading_file & x) {
        throw XProj("Configuration file (monalisa.conf) not found\n");
    }
    notify(vm);

    if (vm.count("help") > 0) {
        cerr << str(format("%s\n") % desc);
        exit(0);
    }

    if (vm.count("version") > 0) {
        cerr << str(format("This is %s version %d.%d (%s)\n") % program_name % major_version % minor_version % minor_version_date);
        cerr << str(format("Written by %s\n") % author);
        exit(0);
    }
    
    if (reference_image_filename.length() == 0) {
        throw XProj("imagefile must be specified");
    }

    if (vm.count("nindivs") > 0) {
        if (Individual::untouchable) {
            if (nindivs < 2) {
                throw XProj("nindivs < 2 makes no sense; this program cannot work with fewer than 2 individuals in its population");
            }
        }
        else {
            if (nindivs < 1) {
                throw XProj("nindivs < 1 makes no sense; this program cannot work with zero individuals in its population");
            }
        }
    }
    
    if (vm.count("nreprod") > 0) {
        if (nreprod > nindivs) {
            throw XProj("nreprod > nindivs makes no sense; can't have more reproducing individuals than individuals in the population");
        }
    }

    if (vm.count("deltargb") > 0) {
        if (Individual::delta_rgb < 0.0) {
            throw XProj("deltargb < 0.0 makes no sense; window size must be positive");
        }
        else if (Individual::delta_rgb > 255.0) {
            throw XProj("deltargb > 255 makes no sense; maximum rgb value is 255");
        }
    }

    if (vm.count("wtmutate") > 0) {
        if (Individual::mutate_weight < 0.0) {
            throw XProj("wtmutate < 0.0 makes no sense; proposal weights must be positive");
        }
    }
}

float updateProgressPlot(unsigned generation, float best_score, CImgDisplay & progdisp) {
    progress.push_back(make_pair(generation, best_score));
    
    // Calculate slope of line connecting last point to previous point
    //  *       ---
    //   \       |     slope is dy/1 = dy
    //    \      |
    //     \     dy
    //      \    |
    //       \   |
    //        * ---
    //  |- 1 -|
    size_t n = progress.size();
    float slope = numeric_limits<float>::max();
    if (n > 6) {
        // average slopes from last point to previous 5 points
        float yultimate = progress[n-1].second;
        float ylag1     = progress[n-2].second;
        float ylag2     = progress[n-3].second;
        float ylag3     = progress[n-4].second;
        float ylag4     = progress[n-5].second;
        float ylag5     = progress[n-6].second;
        float s1 = (yultimate - ylag1)/1.0;
        float s2 = (yultimate - ylag2)/2.0;
        float s3 = (yultimate - ylag3)/3.0;
        float s4 = (yultimate - ylag4)/4.0;
        float s5 = (yultimate - ylag5)/5.0;
        slope = (s1 + s2 + s3 + s4 + s5)/5.0;
    }
    else if (n > 1) {
        float yultimate    = progress[n-1].second;
        float ypenultimate = progress[n-2].second;
        float dy = yultimate - ypenultimate;
        slope = dy;
    }
    
    if (progress.size() > 1) {
        // Define color to use for plot
        const unsigned char line_color[] = { 255, 255, 0 }; // yellow
        
        // Line should be opaque
        float opaque = 1.0;
        
        // Get number of saved data points
        size_t n = progress.size();
        
        // Get dimensions of display
        int w = progdisp.width();
        int h = progdisp.height();
        
        // Create image on which to plot
        CImg<unsigned char> plotimg(w, h, 1, 3, 0);
        
        // Get dimensions of generation and score
        float maxscore = progress[0].second;
        float minscore = progress[n-1].second;
        float firstgen = progress[0].first;
        float lastgen = progress[n-1].first;
        
        // Plot the first line segment
        int x0 = w*(progress[0].first - firstgen)/(lastgen - firstgen);
        int y0 = h - h*(progress[0].second - minscore)/(maxscore - minscore);
        int x1 = w*(progress[1].first - firstgen)/(lastgen - firstgen);;
        int y1 = h - h*(progress[1].second - minscore)/(maxscore - minscore);
        plotimg.fill(0).draw_line(x0, y0, x1, y1, line_color, opaque);
        
        // Plot remaining line segments
        for (unsigned i = 2; i < n; ++i) {
            x0 = x1;
            y0 = y1;
            x1 = w*(progress[i].first - firstgen)/(lastgen - firstgen);;
            y1 = h - h*(progress[i].second - minscore)/(maxscore - minscore);
            plotimg.draw_line(x0, y0, x1, y1, line_color, opaque);
        }
        
        // Show image on display
        plotimg.display(progdisp);
    }
        
    if (progress.size() > 0) {
        ofstream progf("progress.R");

        progf << "gen <- c(" << progress[0].first;
        for (unsigned i = 1; i < progress.size(); ++i)
            progf << "," << progress[i].first;
        progf << ")" << endl;

        progf << "score <- c(" << progress[0].second;
        for (unsigned i = 1; i < progress.size(); ++i)
            progf << "," << progress[i].second;
        progf << ")" << endl;

        progf << "plot(gen, score, type=\"b\", pch=20, lwd=2, col=\"navy\", xlab=\"generation\", ylab=\"best score\")" << endl;

        progf.close();
    }
    return slope;
}

void scoreIndivRange(vector<Individual> & population, vector<pair<float,unsigned> > & scores, const CImg<img_t> & refimg, unsigned first, unsigned last, bool slow) {
    for (unsigned i = first; i < last; ++i) {
        double s = slow ? population[i].calcScore(refimg) :  population[i].calcScoreXY(refimg);
        scores[i].first = s;
        scores[i].second = i;
    }
}

void resetPrevScore(vector<Individual> & population) {
    for (auto it = population.begin(); it != population.end(); ++it) {
        it->resetPrevScore();
    }
}

void scorePopulation(vector<Individual> & population, vector<pair<float,unsigned> > & scores, const CImg<img_t> & refimg, ThreadSchedule & thread_schedule, bool slow) {
    if (nthreads == 1) {
        unsigned i = 0;
        for (auto it = population.begin(); it != population.end(); ++it) {
            double s = slow ? it->calcScore(refimg) : it->calcScoreXY(refimg);
            scores[i].first = s;
            scores[i].second = i;
            ++i;
        }
    }
    else {
        vector<thread> threads;
        for (unsigned i = 0; i < nthreads; i++) {
            threads.push_back(thread(&scoreIndivRange,
                ref(population),
                ref(scores),
                ref(refimg),
                thread_schedule._first[i],
                thread_schedule._last[i],
                slow)
            );
        }

        // The join function causes this loop to pause until the ith thread finishes
        for (unsigned i = 0; i < threads.size(); i++) {
            threads[i].join();
        }
    }
    
    // Sort scores from smallest to largest
    sort(scores.begin(), scores.end());
}

int main(int argc, const char * argv[]) {
    processCommandLineOptions(argc, argv);
    
    ofstream outf(output_filename);
    
    // Set the seed for the random number generator
    lot.setSeed(rnseed);

    // Calculate cumulative mutation probabilites given specified weights
    Individual::calcCumProb();

    // Create image from painting of Mona Lisa
    CImg<img_t> refimg(reference_image_filename.c_str());
        
    // Get width and height of monalisa
    Individual::refw = refimg.width();
    Individual::refh = refimg.height();
    
    // Provide basic information about this run
    if (nthreads > 1) {
        cout << "Using multithreading with " << nthreads << " threads" << endl;
        outf << "Using multithreading with " << nthreads << " threads" << endl;
    }
    else {
        cout << "Using single thread" << endl;
        outf << "Using single thread" << endl;
    }

    if (Individual::untouchable) {
        cout << "Best individual is untouchable" << endl;
        outf << "Best individual is untouchable" << endl;
    }
    else {
        cout << "Best individual can be lost" << endl;
        outf << "Best individual can be lost" << endl;
    }

    cout << "Population size: " << nindivs << endl;
    outf << "Population size: " << nindivs << endl;

    cout << "No. indivs. allowed to reproduce: " << nreprod << endl;
    outf << "No. indivs. allowed to reproduce: " << nreprod << endl;

    ThreadSchedule thread_sched(nthreads, nindivs);
    cout << thread_sched.summary() << endl;
    outf << thread_sched.summary() << endl;

    cout << "reference image width    = " << Individual::refw << endl;
    outf << "reference image width    = " << Individual::refw << endl;

    cout << "reference image height   = " << Individual::refh << endl;
    outf << "reference image height   = " << Individual::refh << endl;

    cout << "reference image spectrum = " << refimg.spectrum() << endl;
    outf << "reference image spectrum = " << refimg.spectrum() << endl;

    CImgDisplay mainwindow;
    CImgDisplay cfwindow;
    CImgDisplay progwindow(Individual::refw, Individual::refh, "Score");

    refimg.display(cfwindow);
        
    // Create a population
    vector<Individual> population(nindivs);
    initializePopulation(population, refimg);
    vector<pair<float,unsigned> > scores(nindivs);
    
    // Score every individual in the population and sort from
    // smallest score (best) to largest score (worst)
    scorePopulation(population, scores, refimg, thread_sched, true);
    resetPrevScore(population);

    // Display and save best individual
    float best_score = scores[0].first;
    float worst_score = scores.rbegin()->first;
    Individual best_indiv = population[scores[0].second];
    best_indiv.display(mainwindow, /*start_message*/true);
    best_indiv.saveToFile("best.jpg", 0);
        
    StopWatch stopwatch;
    stopwatch.start();
    
    // Simulate truncation selection until generation equals stop_at_generation
    // or the user chooses to stop by pressing Q key in any window
    float fitness_slope = numeric_limits<float>::max();
    unsigned generation = 0;
    bool done = false;
    bool running = false;
    while (!done) {
        if (running) {
            generation++;
                    
            // Sanity checks to ensure scores were sorted properly
            assert(scores[0].first <= scores[nreprod-1].first);
            if (nindivs > 2) {
                assert(scores[nreprod-1].first <= scores[nreprod].first);
                assert(scores[nreprod].first <= scores[nindivs-1].first);
            }
            
            // The value of retain[k] will be true if population[k] is one of
            // the nreprod best individuals in the population.
            vector<bool> retain(nindivs, false);
            vector<unsigned> made_the_cut(nreprod);
            for (unsigned i = 0; i < nreprod; ++i) {
                unsigned made_the_cut_index = scores[i].second;
                retain[made_the_cut_index] = true;
                made_the_cut[i] = made_the_cut_index;
            }
            sort(made_the_cut.begin(), made_the_cut.end());
            
            // Replace the individuals that didn't make the cut
            unsigned j = 0;                 // index into made_the_cut
            unsigned k = made_the_cut[j];   // next individual that made the cut
            
            assert(made_the_cut.size() == nreprod);
            unsigned check_number_copied = 0;
            
            for (unsigned child_index = 0; child_index < nindivs; ++child_index) {
                if (child_index < k) {
                    // Choose a parent
                    unsigned pos = lot.randint(0, nreprod - 1);
                    unsigned parent_index = made_the_cut[pos];
    
                    // Copy parent's genome to the child
                    check_number_copied++;
                    population[child_index] = population[parent_index];
                }
                else {
                    if (j < nreprod - 1)
                        k = made_the_cut[++j];
                    else
                        k = nindivs;
                }
            }
            if (check_number_copied != nindivs - nreprod) {
                cerr << "\ngeneration " << generation << ":" << endl;
                
                cerr << "\nretain: ";
                for (unsigned i = 0; i < nindivs; ++i) {
                    cerr << (retain[i] ? "X" : ".");
                }
                cerr << endl;
                
                cerr << "\nmade_the_cut: ";
                for (unsigned i = 0; i < nreprod; ++i) {
                    cerr << made_the_cut[i] << " ";
                }
                cerr << endl;
                throw XProj(format("check_number_copied = %d, nindivs = %d, nreprod = %d, nindivs - nreprod = %d\n ") % check_number_copied % nindivs % nreprod % (nindivs - nreprod));
            }
    
            // Let all individuals undergo mutation, except preserve the best
            // individual if untouchable is true
            unsigned untouchable_index = (Individual::untouchable ? made_the_cut[0] : nindivs);
            for (unsigned i = 0; i < nindivs; ++i) {
                if (i != untouchable_index)
                    population[i].mutate();
            }
            
            scorePopulation(population, scores, refimg, thread_sched, false);
            resetPrevScore(population);
            
            // Check if we have a new all-time best individual
            if (scores[0].first < best_score) {
                best_score = scores[0].first;
                best_indiv = population[scores[0].second];
                if (verbose) {
                    cout << "new best score (" << best_score << ") in generation " << generation << endl;
                    outf << "new best score (" << best_score << ") in generation " << generation << endl;
                }
            }
            
            // Check if we have a new all-time worst score
            if (scores.rbegin()->first > worst_score) {
                worst_score = scores.rbegin()->first;
            }
                    
            if (generation % report_every == 0) {
                fitness_slope = updateProgressPlot(generation, best_score, progwindow);
                string msg = str(format("generation %d\n  best score = %.5f\n  fitness slope = %.5f")
                    % generation
                    % best_score
                    % fitness_slope
                    );
                cout << msg << endl;
                outf << msg << endl;
                
                if (display_all_time_best) {
                    best_indiv.display(mainwindow, /*start_message*/false);
                }
                else {
                    population[0].display(mainwindow, /*start_message*/false);
                }
            }
            
            if (generation % save_every == 0) {
                best_indiv.saveToFile("best.jpg", generation/save_every);
            }
        }

        unsigned keymain = mainwindow.key();
        unsigned keycf   = cfwindow.key();
        unsigned keyplot = progwindow.key();
        if (keymain == cimg::keyQ || keycf == cimg::keyQ || keyplot == cimg::keyQ) {
            best_indiv.saveToFile("best-final.jpg", -1);
            done = true;
        }
        else if (keymain == cimg::keyS || keycf == cimg::keyS || keyplot == cimg::keyS) {
            running = true;
            best_indiv.replaceGenome(::saved_image);
            best_indiv.display(mainwindow, /*start_message*/false);
        }
        
        if (generation == stop_at_generation) {
            best_indiv.saveToFile("best-final.jpg", -1);
            done = true;
        }
    }

    float seconds = stopwatch.stop();

    cout << str(format("Total time: %.3f seconds") % seconds) << endl;
    outf << str(format("Total time: %.3f seconds") % seconds) << endl;

    outf.close();
    
    return 0;
}
