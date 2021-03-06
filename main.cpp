/*
 * PsoPath: Shortest path calculation using Particle Swarm Optimisation
 * Copyright (C) 2020-2021 by Constantine Kyriakopoulos
 * zfox@users.sourceforge.net
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

bool simpleRun()
{
	AdaptiveSystem* pso = new PsoSystem("topology.json", 
			PsoSystem::POPULATION_SIZE, PsoSystem::ITERATIONS);
	auto nodePath = pso->path(0, 19);
	for(int node : nodePath)
		std::cout << node << " ";
	std::cout << std::endl;
	delete pso;

	return nodePath.size() > 0;
}

int main(int argc, char *argv[])
{
	return simpleRun() ? EXIT_SUCCESS : EXIT_FAILURE;
}
