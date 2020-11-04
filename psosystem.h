/*
 * PsoPath: Shortest path calculation using Particle Swarm Optimisation
 * Copyright (C) 2020-2021 by Constantine Kyriakopoulos
 * @version 0.9
 * 
 * @section LICENSE
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef PSOSYSTEM_H
#define PSOSYSTEM_H

#include <iostream>
#include <vector>
#include <set>
#include <exception>
#include <memory>
#include <functional>
#include <random>
#include <limits>
#include <cmath>
#include "adaptivesystem.h"

struct Particle
{
	Particle();
	static long int counter;
	long int id;
	std::vector<std::pair<int, double>> nodePriVec;
	std::vector<double> velocityVec;
	std::vector<double> pBestVec;
	std::vector<double> nBestVec;
	std::vector<int> path;
	std::vector<int> bestPath;
	double fitness;
	double bestFitness;
	bool operator<(const Particle&) const;
	bool operator>(const Particle&) const;
	bool operator==(const Particle&) const;
};

class PsoSystem : public AdaptiveSystem
{
public:
	static const int POPULATION_SIZE = 40;
	static const int ITERATIONS = 150;
	static const int M = 4;
	static const int VELOCITY_LIMIT = 3000;
	static const int VELOCITY_INIT_LIMIT = 10;
	static const int PRIORITY_INIT_LIMIT = 100;
	static constexpr double C1 = 2.05;
	static constexpr double C2 = 2.05;
	static constexpr double FI = C1 + C2;
	PsoSystem(const std::string&, int, int, bool = false);
	virtual ~PsoSystem();
	virtual std::vector<int> path(int, int);
	virtual void clear();

protected:
	virtual void updateFitness(Particle&) const noexcept(false);

private:
	void initParticles(int, int);
	void generatePri(Particle&, int, int);
	void constructPath(Particle&);
	bool isValid(int, int, std::vector<int>&) const;
	bool isNeighbour(int, int) const;
	void updateVelocity(Particle&);
	void updatePosition(Particle&);
	void updatePBest(Particle&);
	void updateNBest(Particle&);
	Particle getPreviousNeighbour(Particle&) const noexcept(false);
	Particle getNextNeighbour(Particle&) const noexcept(false);
	std::vector<int> bestPath() const;
	std::vector<Particle> particles;
	int nodes;
	bool useM;
	int popSize;
	int iterations;
};

#endif // PSOSYSTEM_H
