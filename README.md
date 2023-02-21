#Realtime File System Optimizer Demo
This project is a demo version of a code that was initially created as a realtime file system optimizer. It has been modified and simplified to serve as an educational tool to showcase the use of genetic algorithms.

#Details
The code in this project is composed of a header file genetic.h and genetic.c.

The genetic.h file defines various constants and data structures used in the implementation of the genetic algorithm. It also contains two inline functions, GENOME_TO_SYNC_RATE() and GENOME_TO_READ_AHEAD(), which are responsible for extracting and converting values from the genome.

The biosphere data structure represents the population of the genes, and each gene is stored as an instance of the genome structure. The genes data structure contains the individuals and their fitness, and the biosphere data structure represents the total time and count of the gene pool.

This demo code is intended to serve as a tool to help understand genetic algorithms and how they can be used in optimizing file systems. It is important to note that this demo is not intended for use in a production environment and is only for educational purposes.

#Usage
To use this code, simply download the genetic.h and genetic.c files and include them in your project.

*You will need to redefine the genome structures to accurately represent the variables to optimize.*
