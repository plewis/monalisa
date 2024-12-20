# The Evolution of Mona Lisa

This software provides a graphical illustration of the power of natural selection, especially as it applies to mimicry. It uses a simulation of natural selection in which fitness is determined by similarity to a reference image. Technically it is artificial selection since humans are providing the environment and fitness function, but the underlying process is identical.

It is useful in teaching because it shows that the exquisite [mimicry](https://en.wikipedia.org/wiki/Mimicry) we see in nature can be achieved by random mutation and natural selection.

This work was inspired by an [app created by Roger Johansson](https://rogerjohansson.blog/2008/12/07/genetic-programming-evolution-of-mona-lisa/), but I was encouraged to switch to allowing genomes to be whole images rather than collections of overlapping semi-transparent polygons following a suggestion by Kent Holsinger.

The results are dramatic. Here are three arbitrary images from a single run showing the starting point (random pixel colors) to the final result (300,000 generations) in a population of size 20 in which the fittest 8 are retained as parents for the worst 12 removed each generation.

### Reference image

![Reference image](https://github.com/plewis/monalisa/blob/main/img/mona-lisa-head.jpg?raw=true)

### Starting image

![Starting image](https://github.com/plewis/monalisa/blob/main/img/start.jpg?raw=true)

### Best individual after 300,000 generations

![Image after 300,000 generations](https://github.com/plewis/monalisa/blob/main/img/300k.jpg?raw=true)

### Best individual after 2 million generations

![Final image](https://github.com/plewis/monalisa/blob/main/img/final.jpg?raw=true)
  
The image used as reference (and included in the _test_ directory) is a public domain image of [Leonardo da Vinci's Mona Lisa](https://en.wikipedia.org/wiki/Mona_Lisa#/media/File:Mona_Lisa,_by_Leonardo_da_Vinci,_from_C2RMF_retouched.jpg).

## Control file

The program settings are expected to be in a file named _monalisa.conf_. The default version of this file looks like this:

    rnseed      = 13579
    nindivs     = 20
    nreprod     = 8
    reportevery = 10000
    saveevery   = 10000
    #stopat      = 300000
    deltargb    = 100
    imagefile   = ../img/mona-lisa-head.jpg
    outfile     = output.txt

Comments begin with `#` and a brief description of the options is provided below:

* `rnseed` is the seed for the pseudorandom number generator. Changing this will result in a different simulation but the results will be very similar.
* `nindivs` is the number of haploid "individuals" in the population.
* `nreprod`: the fittest `nreprod` individuals are preserved each generation and serve as parents for the offspring that will replace the `nindivs - nreprod` individuals that are culled each generation. Each offspring is simply a copy of a randomly-chosen parent individual.
* `reportevery` determines how many generations elapse before the status is reported to the console and output file.
* `saveevery` determines how many generations elapse between saved images. That is, the best individual's image is saved every `saveevery` generations.
* `stopat` causes the program to stop at the generation provided. This option is seldom very useful. You probably want to comment this out or remove it entirely and let the program run until you press the Q key in one of the windows.
* `deltargb` determines the maximum distance a mutation to any particular color (red, green, or blue) can be from its starting point. This value must lie between 0 and 255. Larger values cause mutations to be more drastic; smaller values lead to mutations that are less severe but also will slow down the convergence rate.
* `imagefile` specifies the reference image to mimic. This can be any size, but it is best to keep it below 500x500 pixels unless you want to wait a long time for results.
* `outfile` provides the name of the file in which console output is saved.

While this program allows multithreading (you can specify `nthreads` greater than 1 in the _monalisa.conf_ file), you will probably find it does not help because the "fitness" calculation is so cheap that there
is more overhead than processing for the multithreaded version. It may help, however, if the reference image is very large.

## Compiling the program

Note: this software was developed on a MacBook Pro M3 running MacOS 14.6.1 and Xcode 16.2 with [XQuartz 2.8.5](https://www.xquartz.org/) providing the X11 windowing system. It should work on other systems but it has not been tested on other systems.

You will need to install the following libraries in order to compile it:

* [CImg](https://cimg.eu/index.html)
* [Boost](https://www.boost.org) [program_options](https://www.boost.org/doc/libs/1_87_0/doc/html/program_options.html)
* [Boost](https://www.boost.org) [format](https://www.boost.org/doc/libs/1_87_0/libs/format/)
* [Boost](https://www.boost.org) [random](https://www.boost.org/doc/libs/1_87_0/doc/html/boost_random.html)

Note that the Boost program options library must be compiled (it is not a header-only library like many Boost libraries).

You will also need to install X11 for the graphical display if it does not come with your operating system.
A makefile is provided that can be adjusted for your system.


