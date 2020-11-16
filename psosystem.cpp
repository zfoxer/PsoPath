/*
 * PsoPath: Shortest path calculation using Particle Swarm Optimisation
 * Copyright (C) 2020 by Constantine Kyriakopoulos
 * @version 1.0.1
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

#include "psosystem.h"

PsoSystem::PsoSystem(const std::string& filename, int popSize, int iterations, bool useM)
{
	initTopo(filename);
	this->useM = useM;
	
	std::set<int> uniqueNodes;
	for(auto& edge : edges)
		uniqueNodes.insert({edge.edgeStart, edge.edgeEnd});
	nodes = uniqueNodes.size();

	if(popSize > 0 && iterations > 0)
	{
		this->popSize = popSize;
		this->iterations = iterations;
	}
	else
	{
		this->popSize = POPULATION_SIZE;
		this->iterations = ITERATIONS;
	}

	std::random_device rd;
	gen = std::mt19937_64(rd());
}

PsoSystem::~PsoSystem() { }

std::vector<int> PsoSystem::path(int src, int dest)
{
	if(nodes <= 1)
		return std::vector<int>();

	initParticles(src, dest);

	for(auto& particle : particles)
	{
		constructPath(particle);
		updateFitness(particle);
		updatePBest(particle);
		updateNBest(particle);                   
	}
	
	int i = 0;
	do
	{
		for(auto& particle : particles)
		{
			updateVelocity(particle);
			updatePosition(particle);
			constructPath(particle);
			updateFitness(particle);
			updatePBest(particle);
			updateNBest(particle);                    
		}       
	}
	while(i++ < iterations);

	return bestPath();
}

void PsoSystem::initParticles(int src, int dest)
{
	particles.clear();
	std::uniform_int_distribution<> distro(-VELOCITY_INIT_LIMIT, VELOCITY_INIT_LIMIT);

	int i = 0;
	while(i++ < popSize)
	{
		Particle tempParticle;
		generatePri(tempParticle, src, dest);
		for(int j = 0; j < nodes; ++j)
			tempParticle.velocityVec.push_back(distro(gen));
		for(int k = 0; k < tempParticle.nodePriVec.size(); ++k)
			tempParticle.pBestVec.push_back(tempParticle.nodePriVec.at(k).second);
		tempParticle.fitness = 0;
		tempParticle.bestFitness = 0;
		
		particles.push_back(tempParticle);
	}

	for(auto& particle : particles)
		updateNBest(particle);
}

void PsoSystem::generatePri(Particle& particle, int src, int dest)
{
	particle.nodePriVec.clear();
	std::uniform_int_distribution<> distro(-PRIORITY_INIT_LIMIT, PRIORITY_INIT_LIMIT);
		
	particle.nodePriVec.push_back({src, distro(gen)});
	for(int node = 0; node < nodes; ++node)
		if(node != src && node != dest)
			particle.nodePriVec.push_back({node, distro(gen)});
	particle.nodePriVec.push_back({dest, distro(gen)});
}

void PsoSystem::constructPath(Particle& particle)
{
	std::vector<int> path;
	std::vector<std::pair<int, double>> npCopy(particle.nodePriVec);
	path.push_back((*npCopy.cbegin()).first);
	int highestPri = (*npCopy.cbegin()).second;
	int highestPriNode = (*npCopy.cbegin()).first;
	(*npCopy.begin()).second = std::numeric_limits<int>::min();

	while(true)
	{
		auto it = npCopy.begin();
		auto tempIt = npCopy.begin();
		while(it != npCopy.end())
		{
			auto ndTemp = (*it).first;
			auto priTemp = (*it).second;
			if(priTemp > highestPri)
			{
				highestPri = priTemp;
				highestPriNode = ndTemp;
				tempIt = it;
			}
			++it;
		}
		
		if(useM)
		{
			if(isNeighbour(highestPriNode, *path.crbegin()) && highestPriNode - *path.crbegin() > -M)
				path.push_back(highestPriNode);
		}
		else
			if(isNeighbour(highestPriNode, *path.crbegin()))
				path.push_back(highestPriNode);

		if(highestPriNode == (*npCopy.crbegin()).first)
		{
			if(isValid((*npCopy.cbegin()).first, (*npCopy.crbegin()).first, path))
				particle.path = path;
			break;
		}
		
		(*tempIt).second = std::numeric_limits<int>::min();
		highestPri = (*npCopy.cbegin()).second;
		highestPriNode = (*npCopy.cbegin()).first;
	}
}

void PsoSystem::updateFitness(Particle& particle) const noexcept(false)
{
	if(!particle.path.size() || particle.path.size() == 1)
	{    
		particle.fitness = 0;
		return;
	}

	double cost = 0;
	std::vector<int>::const_iterator it = particle.path.cbegin();
	while(it++ != particle.path.cend() - 1)
	{
		int edgeStart = *(it - 1);
		int edgeEnd = *it;

		auto tempIt = std::find_if(edges.cbegin(), edges.cend(), \
				[edgeStart, edgeEnd](Edge edge)
				{
					return edge.edgeStart == edgeStart && edge.edgeEnd == edgeEnd;
				});
		if(tempIt == edges.cend())
			throw std::invalid_argument("\nPsoSystem::updateFitness... Edge not found");

		cost += (*tempIt).weight;
	}

	particle.fitness = 1 / cost;
}

void PsoSystem::updateVelocity(Particle& particle)
{
	std::uniform_real_distribution<> distro(0, 1);
	
	for(int i = 0; i < nodes; ++i)
	{
		double r1 = distro(gen);
		double r2 = distro(gen);
		double tempVelo = (C1 * r1 * (particle.pBestVec.at(i) - particle.nodePriVec.at(i).second) + C2 * r2 * 
				(particle.nBestVec.at(i) - particle.nodePriVec.at(i).second));
	
		tempVelo += particle.velocityVec.at(i);
		double xi = 2 / fabs((2 - FI - sqrt(FI * FI - 4 * FI)));
	
		double velocity = xi * tempVelo;
		if(velocity > VELOCITY_LIMIT)
			velocity = VELOCITY_LIMIT;
		if(velocity < -VELOCITY_LIMIT)
			velocity = -VELOCITY_LIMIT;

		particle.velocityVec[i] = velocity;
	}
}

void PsoSystem::updatePosition(Particle& particle)
{
	for(int i = 0; i < nodes; ++i)
		particle.nodePriVec[i].second = particle.nodePriVec[i].second + particle.velocityVec.at(i);
}

void PsoSystem::updatePBest(Particle& particle)
{
	if(particle.fitness > particle.bestFitness)
	{
		particle.bestFitness = particle.fitness;
		particle.bestPath = particle.path;
		for(int i = 0; i < nodes; ++i)
			particle.pBestVec[i] = particle.nodePriVec.at(i).second; 
	}
}
	
void PsoSystem::updateNBest(Particle& particle)
{
	auto nextNeighbour = getNextNeighbour(particle);
	auto previousNeighbour = getPreviousNeighbour(particle);
	particle.nBestVec = nextNeighbour.bestFitness > previousNeighbour.bestFitness 
			? nextNeighbour.pBestVec : previousNeighbour.pBestVec;
}

std::vector<int> PsoSystem::bestPath() const
{
	std::vector<int> bestPath;
	double tempValue = std::numeric_limits<double>::min();
	for(auto& particle : particles)
		if(particle.bestFitness > tempValue)
		{
			tempValue = particle.bestFitness;
			bestPath = particle.bestPath;
		}

	return bestPath;
}

Particle PsoSystem::getPreviousNeighbour(Particle& particle) const noexcept(false)
{
	if(particles.size() <= 1)
		throw std::out_of_range("\nPsoSystem::getPreviousNeighbour... No neighbour");

	auto it = std::find_if(particles.cbegin(), particles.cend(), 
			[particle](Particle prt)
			{
				return particle.id == prt.id;
			});

	if(it == particles.cend())
		throw std::invalid_argument("\nPsoSystem::getPreviousNeighbour... Particle not found");

	return (it == particles.cbegin()) ? *particles.crbegin() : *(it - 1);
}
	
Particle PsoSystem::getNextNeighbour(Particle& particle) const noexcept(false)
{
	if(particles.size() <= 1)
		throw std::out_of_range("\nPsoSystem::getNextNeighbour... No neighbour");

	auto it = std::find_if(particles.cbegin(), particles.cend(), 
			[particle](Particle prt)
			{
				return particle.id == prt.id;
			});

	if(it == particles.cend())
		throw std::invalid_argument("\nPsoSystem::getNextNeighbour... Particle not found");

	return (it == (particles.cend() - 1)) ? *particles.cbegin() : *(it + 1);
}

bool PsoSystem::isValid(int src, int dest, std::vector<int>& path) const
{
	return (*path.begin() != src || *path.rbegin() != dest) ? false : true;
}

bool PsoSystem::isNeighbour(int a, int b) const
{
	for(auto& edge : edges)
		if((edge.edgeStart == a && edge.edgeEnd == b) || (edge.edgeStart == b && edge.edgeEnd == a))
			return true;

	return false;
}

void PsoSystem::clear()
{
	edges.clear();
	particles.clear();
	useM = nodes = 0;
}

Particle::Particle()
{
	id = ++counter;
	fitness = bestFitness = 0;
}

bool Particle::operator<(const Particle& rhs) const
{
	return id < rhs.id;
}

bool Particle::operator>(const Particle& rhs) const
{
	return id > rhs.id;
}

bool Particle::operator==(const Particle& rhs) const
{
	return id == rhs.id;
}

long int Particle::counter = 0;
