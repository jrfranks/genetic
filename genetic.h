#ifndef _GENETIC_H
#define _GENETIC_H

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
 * @brief This file contains definitions and functions for algorithms in
 * genetic.c.
 *
 * The genetic algorithm is a heuristic optimization algorithm inspired by the
 * process of natural selection. It uses mechanisms such as reproduction,
 * mutation, and selection to find the optimal solution to a problem. This file
 * contains constants and data structures used by the genetic algorithm, as well
 * as two conversion functions that convert genome values to sync rate and read
 * ahead values, respectively. The genome is a data structure that represents
 * the genetic information of an individual in the population. The genes data
 * structure represents the entire population of individuals, and the biosphere
 * data structure contains information about the execution of the genetic
 * algorithm.
 */

#define POPULATION_SIZE 32

#define GENOME_SYNC_BITS 10
#define GENOME_RA_BITS 2
#define GENOME_BITS_USED (GENOME_SYNC_BITS + GENOME_RA_BITS)

#define GENOME_SYNC_MASK ((1 << GENOME_SYNC_BITS) - 1)
#define GENOME_RA_MASK (((1 << GENOME_RA_BITS) - 1) << GENOME_SYNC_BITS)
#define GENOME_BITS_MASK ((1 << GENOME_BITS_USED) - 1)

#define INIT_GENOME(s, r)                                                      \
	(((s)&GENOME_SYNC_MASK) | ((r << GENOME_SYNC_BITS) & GENOME_RA_MASK))

// Structure definitions

/*
 * Genome and the current fitness of an individual.
 */
typedef struct genome {
	uint64_t		fitness;	// Current fitness score
	uint16_t		genome;		// of this genome.
} genome_t;

/*
 * This structure tracks all the individuals in the gene pool
 * along with enough information to give each intividual a
 * chance at reporduction.
 */
typedef struct genes {
	int16_t		current_individual;		// Index of individual being tested
	int16_t		cycles_remaining;		// Number of test measurements remaining.
	int16_t		roulette_slot;			// 
	genome_t	individuals[POPULATION_SIZE];
} genes_t;

/*
 * This structure tracks the trial progress of the enclosed gene pool.
 */
typedef struct biosphere {
	uint64_t	total_time;
	uint32_t	total_count;
	genes_t		gene_pool;
} biosphere_t;


/*
 * Root structure for optimization
 *
 * For this implementation we are optimizing sync_rate and read ahead.
 *
 * All the trial data is stored in the biosphere structure.
 */

typedef struct online {
	char *		name;
	biosphere_t	biosphere;
	uint32_t	sync_rate;
	uint32_t	read_ahead;
} online_t;

// Funtion prototypes
void	initialize_genetic_pools(online_t *fs);
void	compute_fitness(online_t *fs);
void	crossover_and_mutate(genome_t *child, const genome_t *mate1,
		const genome_t *mate2, int32_t mutation_rate);
genome_t *select_mate(genes_t *gp, uint32_t total_fitness);
void	evolve(online_t *fs);
void	optimize(online_t *fs, uint64_t response_time);
void	stats(online_t *fs);

#endif /* !_GENETIC_H */
