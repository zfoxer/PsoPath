/*
 * PsoPath: Shortest path calculation using Particle Swarm Optimisation
 * Copyright (C) 2020 by Constantine Kyriakopoulos
 * @version 1.0.2
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

/**
 * Class constructor for the PSO system.
 *  
 * @param filename File containing the JSON representation of the topology
 * @param popSize Population size for the PSO system
 * @param iterations Number of iterations while constructing path
 * @param useM If the M parameter is considered, looking back only M places while traversing
 */
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

/**
 * Class constructor for the PSO system starting with empty topology.
 *  
 * @param popSize Population size for the PSO system
 * @param iterations Number of iterations while constructing path
 * @param useM The M parameter is considered, looking back only M places while traversing
 */
PsoSystem::PsoSystem(int popSize, int iterations, bool useM)
{
	this->useM = useM;
	
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

/**
 *  Destructor.
 */
PsoSystem::~PsoSystem() { }

/**
 * Finds the path between src and dest.
 *
 * @param src Originating node
 * @param dest Destination node of the request
 * @return std::vector<int> Node container representing the path
 */
std::vector<int> PsoSystem::path(int src, int dest)
{
	if(nodes <= 1)
		return std::vector<int>();

	// Initialise all particles in their vector (velocity & position)
	initParticles(src, dest);

	// For each particle, construct path and store its fitness value
	// Iterate and update pBest and nBest in each particle
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
		// Reposition particles
		for(auto& particle : particles)
		{
			// Move particle into the search space
			updateVelocity(particle);
			updatePosition(particle);
			// Construct path and evaluate fitness of each particle
			constructPath(particle);
			updateFitness(particle);
			// Update pBest & nBest vectors
			updatePBest(particle);
			updateNBest(particle);                    
		}       
	}
	while(i++ < iterations);

	return bestPath();
}

/**
 * Initialises the particle population.
 *
 * @param src Originating node
 * @param dest Destination node of the request
 */
void PsoSystem::initParticles(int src, int dest)
{
	particles.clear();
	std::uniform_int_distribution<> distro(-VELOCITY_INIT_LIMIT, VELOCITY_INIT_LIMIT);

	// Priorities and velocities
	int i = 0;
	while(i++ < popSize)
	{
		// Particle creation
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

/**
 * Generates the particle priority vector.
 *
 * @param particle The particle to generate its node priority vector
 * @param src Originating node
 * @param dest Destination node of the request
 */
void PsoSystem::generatePri(Particle& particle, int src, int dest)
{
	particle.nodePriVec.clear();
	std::uniform_int_distribution<> distro(-PRIORITY_INIT_LIMIT, PRIORITY_INIT_LIMIT);
		
	// Topology nodes mapped to priorities
	particle.nodePriVec.push_back({src, distro(gen)});
	for(int node = 0; node < nodes; ++node)
		if(node != src && node != dest)
			particle.nodePriVec.push_back({node, distro(gen)});
	particle.nodePriVec.push_back({dest, distro(gen)});
}

/**
 * Constructs the internal path of the given particle.
 *
 * @param particle This particle's path that will be constructed
 */
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
		// Find the node with the highest priority
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
		
		// It should be a neighbour to the previously inserted node 
		// and less from the previous (node number subtraction) no more than M
		if(useM)
		{
			// Explicit braces avoid dangling 'else'
			if(isNeighbour(highestPriNode, *path.crbegin()) && highestPriNode - *path.crbegin() > -M)
				path.push_back(highestPriNode);
		}
		else
			if(isNeighbour(highestPriNode, *path.crbegin()))
				path.push_back(highestPriNode);

		// Check the last one
		if(highestPriNode == (*npCopy.crbegin()).first)
		{
			if(isValid((*npCopy.cbegin()).first, (*npCopy.crbegin()).first, path))
				particle.path = path;
			break;
		}
		
		// The chosen node will now have the lowest possible priority, won't be chosen again
		(*tempIt).second = std::numeric_limits<int>::min();
		highestPri = (*npCopy.cbegin()).second;
		highestPriNode = (*npCopy.cbegin()).first;
	}
}

/**
 * Updates particle's fitness. In this case, the path cost.
 *
 * @param particle The particle to update its fitness
 */
void PsoSystem::updateFitness(Particle& particle)
{
	// It should include a valid path
	if(!particle.path.size() || particle.path.size() == 1)
	{    
		particle.fitness = 0;
		return;
	}

	// For every path's edge, find its cost
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

/**
 * Updates particle's velocity.
 *
 * @param particle The particle to update its velocity
 */
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
		// Clamp velocity
		if(velocity > VELOCITY_LIMIT)
			velocity = VELOCITY_LIMIT;
		if(velocity < -VELOCITY_LIMIT)
			velocity = -VELOCITY_LIMIT;

		particle.velocityVec[i] = velocity;
	}
}

/**
 * Updates particle's position.
 *
 * @param particle The particle to update its position
 */
void PsoSystem::updatePosition(Particle& particle)
{
	for(int i = 0; i < nodes; ++i)
		particle.nodePriVec[i].second = particle.nodePriVec[i].second + particle.velocityVec.at(i);
}

/**
 * Updates particle's pBest value.
 *
 * @param particle The particle to update its pBest value, i.e., the best position it had so far
 */
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
	
/**
 * Updates particle's nBest value.
 *
 * @param particle Particle to update its nBest value: the best pBest value of its neighbours
 */
void PsoSystem::updateNBest(Particle& particle)
{
	auto nextNeighbour = getNextNeighbour(particle);
	auto previousNeighbour = getPreviousNeighbour(particle);
	particle.nBestVec = nextNeighbour.bestFitness > previousNeighbour.bestFitness 
			? nextNeighbour.pBestVec : previousNeighbour.pBestVec;
}

/**
 * Returns the best path checking all particles.
 *
 * @return std::vector<int> The node container representing the path
 */
std::vector<int> PsoSystem::bestPath() const
{
	// The path with the highest fitness value
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

/**
 * Returns the previous neighbour of a ring particle topology.
 *
 * @param particle The particle to return its previous neighbour
 * @return Particle The previous neighbour
 */
Particle PsoSystem::getPreviousNeighbour(Particle& particle) const
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
	
/**
 * Returns the next neighbour of a ring particle topology.
 *
 * @param particle The particle to return its next neighbour
 * @return Particle The next neighbour
 */
Particle PsoSystem::getNextNeighbour(Particle& particle) const
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

/**
 * Checks the validity of the given path.
 *
 * @param src Path's source node
 * @param dest Path's destination node
 * @param path The node path to validate
 * @return bool Validity indicator
 */
bool PsoSystem::isValid(int src, int dest, std::vector<int>& path) const
{
	// Ensure path's src and dest are correct
	return (*path.begin() != src || *path.rbegin() != dest) ? false : true;
}

/**
 * Checks if the given nodes are neighbours.
 *
 * @param a The first node
 * @param b The second node
 * @return bool True if nodes are neighbours
 */
bool PsoSystem::isNeighbour(int a, int b) const
{
	for(auto& edge : edges)
		if((edge.edgeStart == a && edge.edgeEnd == b) || (edge.edgeStart == b && edge.edgeEnd == a))
			return true;

	return false;
}

/**
 * Clears the system.
 */
void PsoSystem::clear()
{
	edges.clear();
	particles.clear();
	useM = nodes = 0;
}

/**
 * Inserts an edge.
 * 
 * @param src Source node
 * @param dest Destination node
 * @param weight Weight for the edge
 */
void PsoSystem::insertEdge(int src, int dest, double weight)
{
	AdaptiveSystem::insertEdge(src, dest, weight);
	std::set<int> uniqueNodes;
	for(auto& edge : edges)
		uniqueNodes.insert({edge.edgeStart, edge.edgeEnd});
	nodes = uniqueNodes.size();
}

/**
 * Constructor for particles.
 */
Particle::Particle()
{
	id = ++counter;
	fitness = bestFitness = 0;
}

/**
 * Comparison of current instance with the rhs, based on ids.
 *
 * @param rhs The right-hand side object
 * @return bool The indication of current id being less than rhs'
 */
bool Particle::operator<(const Particle& rhs) const
{
	return id < rhs.id;
}

/**
 * Comparison of current instance with the rhs, based on ids.
 *
 * @param rhs The right-hand side object
 * @return bool The indication of current id being greater than rhs'
 */
bool Particle::operator>(const Particle& rhs) const
{
	return id > rhs.id;
}

/**
 * Comparison of current instance with the rhs for equality, based on ids.
 *
 * @param rhs The right-hand side object
 * @return bool The indication of equality
 */
bool Particle::operator==(const Particle& rhs) const
{
	return id == rhs.id;
}

/**
 * Used for producing particle IDs.
 */
long int Particle::counter = 0;
