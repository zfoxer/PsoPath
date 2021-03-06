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

#ifndef ADAPTIVESYSTEM_H
#define ADAPTIVESYSTEM_H

#include <functional>
#include <vector>
#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

using boost::property_tree::ptree;

class AdaptiveSystem
{
public:
	struct Edge
	{
		Edge();
		int edgeStart;
		int edgeEnd;
		double weight;
		long int id;
		bool operator<(const AdaptiveSystem::Edge&) const;
		bool operator>(const AdaptiveSystem::Edge&) const;
		bool operator==(const AdaptiveSystem::Edge&) const;
	};

	AdaptiveSystem();
	virtual ~AdaptiveSystem();
	virtual std::vector<int> path(int, int) = 0;
	virtual void insertEdge(int, int, double) noexcept(false);
	virtual void clear() = 0;

protected:
	virtual void initTopo(const std::string&);
	std::vector<AdaptiveSystem::Edge> edges;

private:
	static int edgeIdCnt;
};

#endif // ADAPTIVESYSTEM_H
