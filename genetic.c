/******************************************************************************
 * Copyright (c) 2022 Sveltesoft. All rights reserved.
 *
 * This software is the confidential and proprietary information of SvelteSoft.
 * You shall not disclose such confidential information and shall use it only in
 * accordance with the terms of the license agreement you entered into with
 * SvelteSoft.
 *****************************************************************************/

/**
 * @file genetic.h
 *
 * @brief This file contains definitions and functions related to genetic
 * algorithms. The genetic algorithm is a heuristic optimization algorithm
 * inspired by the process of natural selection. It uses mechanisms such as
 * reproduction, mutation, and selection to find the optimal solution to a
 * problem.
 *
 * This file contains the source code for a genetic algorithm designed to
 * optimize a file system in realtime. The functions in the file perform various
 * tasks related to the genetic algorithm, including initializing genetic pools,
 * computing fitness, optimizing file systems, and displaying statistics. The
 * file includes detailed function comments in Doxygen style to describe the
 * purpose and functionality of each function.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "genetic.h"

// Initial file system configuration values.
int32_t max_sync_rate = 100;
int32_t min_sync_rate = 5;
int32_t max_read_ahead = 32;

int32_t ops_per_trial = 200;
uint32_t mutation_rate = 20;

/**
 * @brief Converts a genome value to a sync rate value within the valid range.
 * The function takes a 16-bit genome value as input and extracts the lower 10
 * bits representing the sync rate. The function then checks if the extracted
 * value is within the range of valid sync rates, which are specified by the
 * global variables max_sync_rate and min_sync_rate. If the value is out of
 * range, it is replaced by the closest valid value. Finally, the function
 * returns the sync rate value.
 *
 * @param genome A 16-bit unsigned short value representing the genome to be
 * converted.
 * @return A long integer representing the sync rate value within the valid
 * range.
 */

static inline long genome_to_sync_rate(uint16_t genome)
{
	genome &= GENOME_SYNC_MASK;
	if (genome > max_sync_rate)
		genome = max_sync_rate;
	else if (genome < min_sync_rate)
		genome = min_sync_rate;
	return genome;
}

/**
 * @brief Converts the read ahead value stored in the genome to an integer value
 * This function takes a 16-bit genome and extracts the read ahead value, which
 * is stored in bits 10-11 of the genome, and returns it as an integer value. If
 * the extracted value is greater than the specified maximum read ahead value,
 * then the maximum read ahead value is returned instead.
 *
 * @param genome A 16-bit unsigned short integer representing a genome
 * @param max_read_ahead An integer representing the maximum allowed read ahead
 * value
 * @return An integer representing the read ahead value stored in the genome
 */

static inline long genome_to_read_ahead(uint16_t genome)
{
	genome = ((genome)&GENOME_RA_MASK) >> GENOME_SYNC_BITS;
	return ((genome > max_read_ahead) ? max_read_ahead : genome);
}

/**
 * @brief Initializes the genetic pools for a specified online file system.
 * This function initializes the genetic pools for the given online file system
 * by setting up the biosphere, gene pool, and individual genomes. It uses the
 * clock to seed the random number generator, clears the biosphere using bzero,
 * and initializes the gene pool cycles remaining to the number of operations
 * per trial. It then generates a population of individuals, each with a
 * randomly generated genome and a fitness value of 0. The function sets the
 * stop rate for the file system, sync rate, and read ahead parameters based on
 * the genome of the first individual in the gene pool.
 *
 * @param fs A pointer to the online file system struct to initialize genetic
 * pools for.
 * @return void
 */

void initialize_genetic_pools(online_t *fs)
{
	genes_t *gp = &fs->biosphere.gene_pool;
	genome_t *ip;

	srand(time(NULL)); // Seed random number generator.
	memset(&fs->biosphere, 0,
	    sizeof(fs->biosphere)); // Clear the biosphere struct.
	gp->cycles_remaining = ops_per_trial;

	for (int32_t i = 0; i < POPULATION_SIZE; i++) {
		ip = &gp->individuals[i];
		ip->genome = (uint16_t)(rand() & GENOME_BITS_MASK);
		ip->fitness = 0;
	}

	ip = &gp->individuals[0];
	fs->sync_rate = genome_to_sync_rate(ip->genome);
	fs->read_ahead = genome_to_read_ahead(ip->genome);
}

/**
 * @brief Computes the fitness of the current individual in the gene pool of the
 * specified online file system using the response time data collected during
 * the previous trial. This function takes a pointer to an "online" struct
 * representing an online file system and computes the fitness of the current
 * individual in the gene pool. The function updates the fitness value of the
 * current individual by dividing the total time by the total count of
 * operations during the previous trial. It then resets the total time and count
 * of the biosphere to zero to prepare for the next trial. The fitness of the
 * current individual is used in selection during the genetic algorithm
 * iteration to determine which individuals will be chosen to reproduce and
 * generate the next generation.
 *
 * @param fs Pointer to the online file system to compute the fitness for.
 * @return void
 */

void compute_fitness(online_t *fs)
{
	genes_t *gp = &fs->biosphere.gene_pool;
	genome_t *ip = &gp->individuals[gp->current_individual];
	uint32_t cnt = fs->biosphere.total_count;

	ip->fitness = 0;
	ip->fitness += fs->biosphere.total_time / cnt;
	fs->biosphere.total_time = 0;
	fs->biosphere.total_count = 0;
}

/**
 * @brief Performs crossover and mutation on two parent genomes to produce a
 * child genome This function takes two parent genomes, mate1 and mate2, and
 * performs crossover and mutation on them to create a child genome. Crossover
 * occurs at two points, one for sync rate and one for read ahead, and is
 * determined randomly for each point. The resulting child genome has the sync
 * rate from one parent and the read ahead from the other. After crossover,
 * there is a chance for mutation to occur, which happens randomly with a
 * probability of mutation_rate. If mutation occurs, a random bit of the child
 * genome is flipped, and the mutation count is incremented.
 *
 * @param mate1 Pointer to the first parent genome
 * @param mate2 Pointer to the second parent genome
 * @param child Pointer to the child genome that will be generated
 * @return void
 */

void crossover_and_mutate(genome_t *child, const genome_t *mate1,
    const genome_t *mate2, int32_t mutation_rate)
{
	uint16_t g1 = mate1->genome;
	uint16_t g2 = mate2->genome;
	uint16_t crossover1, crossover2;
	uint16_t mask;
	uint16_t gc;

	/* Crossover Sync Rate */
	crossover1 = rand() % GENOME_SYNC_BITS;
	mask = (1 << crossover1) - 1;

	/* Crossover Read Ahead */
	crossover2 = (rand() % GENOME_RA_BITS);
	mask |= ((1 << crossover2) - 1) << GENOME_SYNC_BITS;

	gc = (g1 & mask) | (g2 & ~mask);

	/* Now mutate */
	if ((rand() % mutation_rate) == 0) {
		gc ^= 1 << (rand() % GENOME_BITS_USED);
	}

	child->genome = gc;
}

/**
 * @brief Selects a mate from the gene pool using a roulette wheel selection
 * method. The function uses a roulette wheel selection method. It chooses the
 * individual with the least execution time as the most probable winner. The selection
 * is based on the total fitness of the gene pool and a random spin.
 *
 * @param gp Pointer to the gene pool.
 * @param total_fitness Total fitness of the gene pool.
 * @return A pointer to the selected mate's genome.
 *
 * @note This function modifies the roulette_slot member of the gene pool.
 */

genome_t *select_mate(genes_t *gp, uint32_t total_fitness)
{
	int32_t slot = gp->roulette_slot;
	int32_t spin = lrand48() % total_fitness;

	/*
	 * Use roulette wheel but use the least contributor to the
	 * total as the most probable winner.  Least response time wins!
	 */
	for (;;) {
		if ((spin -= total_fitness - gp->individuals[slot].fitness) < 0)
			break;

		/* Still spinning... */
		slot = (slot + 1) % POPULATION_SIZE;
	}
	gp->roulette_slot = (slot + 1) % POPULATION_SIZE;
	return &gp->individuals[slot];
}

/**
 * @brief Performs a genetic algorithm iteration on the gene pool of the given
 * online file system. This function takes a pointer to an "online" struct
 * representing an online file system and evolves its gene pool using a genetic
 * algorithm. It performs selection, crossover, and mutation to generate a new
 * generation of individuals, which replaces the current generation in the gene
 * pool.
 *
 * @param fs Pointer to an "online" struct representing an online file system to
 * optimize.
 * @return None
 *
 * @note This function modifies the gene pool of the online file system passed
 * as argument.
 */

void evolve(online_t *fs)
{
	genes_t new_gene_pool;
	genes_t *gp = &fs->biosphere.gene_pool;

	// Sum up fitness values
	uint32_t total_fitness = 0;
	for (int32_t i = 0; i < POPULATION_SIZE; i++) {
		total_fitness += gp->individuals[i].fitness;
	}

	// Create a new gene pool
	new_gene_pool.cycles_remaining = gp->cycles_remaining;
	for (int32_t i = 0; i < POPULATION_SIZE; i++) {
		genome_t *mate1 = select_mate(gp, total_fitness);
		genome_t *mate2 = select_mate(gp, total_fitness);
		crossover_and_mutate(mate1, mate2, &new_gene_pool.individuals[i],
				mutation_rate);
	}

	// Replace old gene pool with new one
	*gp = new_gene_pool;
}

/**
 * @brief Performs genetic algorithm optimization for a specified online
 * filesystem. This function evolves the population of the given online
 * filesystem using a genetic algorithm to optimize its parameters. The response
 * time for the current trial is passed as a parameter. The function updates the
 * total time and count for the biosphere, and calculates the fitness for the
 * current individual. If all trials for the current individual have completed,
 * the function evolves the population and resets the fitness counters for the
 * next set of trials. The function then loads the parameters for the next trial
 * from the genome of the current individual.
 * @param fs Pointer to the online filesystem struct to optimize.
 *
 * @param response_time The response time for the current trial.
 * @return void
 *
 * @note The function modifies the population of the online filesystem.
 */

void optimize(online_t *fs, uint64_t response_time)
{
	genes_t *active = &fs->biosphere.gene_pool;
	genome_t *ip;

	// Update trial statistics
	fs->biosphere.total_time += response_time;
	fs->biosphere.total_count++;

	// Check if all trials for current individual are completed
	if (--active->cycles_remaining <= 0) {
		// Calculate fitness and move to next individual
		active->cycles_remaining = ops_per_trial;
		compute_fitness(fs);
		active->current_individual =
		    (active->current_individual + 1) % POPULATION_SIZE;

		// Evolve if all individuals have been tried
		if (active->current_individual == 0) {
			// All trials are over...  Time to evolve.
			evolve(fs);

			// Reset fitness counters for new trials
			for (int32_t i = 0; i < POPULATION_SIZE; i++) {
				active->individuals[i].fitness = 0;
			}
		}

		// Load new individual's parameters for next trial
		ip = &active->individuals[active->current_individual];
		fs->sync_rate = genome_to_sync_rate(ip->genome);
		fs->read_ahead = genome_to_read_ahead(ip->genome);
	}
}

/**
 * @brief Displays statistics for genetic algorithm runs on a specified
 * filesystem or all filesystems. This function displays various statistics for
 * genetic algorithm runs on the specified filesystem or on all filesystems if
 * name is NULL. The statistics include sync rate, read ahead, gene pool
 * information, and genome information for each individual in the population.
 * The average sync rate and read ahead are also calculated and displayed.
 *
 * @param name The name of the filesystem to display statistics for, or NULL for
 * all filesystems.
 * @return 0 on success.
 */

void stats(online_t *fs)
{
	genes_t *gp;
	genome_t *ip;
	uint32_t avg_sync = 0, avg_ra = 0;
	int32_t i;

	printf("%s: sync_rate=%d, read_ahead=%d\n", fs->name,
	    fs->sync_rate, fs->read_ahead);

	gp = &fs->biosphere.gene_pool;
	printf("  cur_ind=%d  cycles_remaining=%d\n",
	    gp->current_individual, gp->cycles_remaining);

	avg_sync = 0;
	avg_ra = 0;
	for (i = 0; i < POPULATION_SIZE; i++) {
		ip = &gp->individuals[i];
		avg_sync += genome_to_sync_rate(ip->genome);
		avg_ra += genome_to_read_ahead(ip->genome);
		printf(
		    "    individual %02d: g=0x03%x, syn=%03ld ra=%ld\n",
		    i, ip->genome, genome_to_sync_rate(ip->genome),
		    genome_to_read_ahead(ip->genome));
	}
	printf("    summary: avg_syn=%d avg_ra=%d\n",
	    avg_sync / POPULATION_SIZE, avg_ra / POPULATION_SIZE);
}
