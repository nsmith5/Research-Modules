#
# Copyright 2016 Nathan Smith
#
# This file is part of Research Modules.
#
# Research Modules is free software: you can redistribute it
# and/or modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# Research Modules is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with Model A Simulation. If not, see http://www.gnu.org/licenses/.
#

module PhaseDiagrams

using Functional, PointGroups
export ModeAnsatz, A₁projection, mode

### Mode expansions ###

type ModeAnsatz <: Function
	q::Vector
	F::Function
	function ModeAnsatz(q::Vector)
		f = (x::Vector) -> cos(q⋅x)
		new(q, f)
	end
end

function mode(f::ModeAnsatz)
	return f.q
end

function (f::ModeAnsatz)(x::Vector)
	return f.F(x)
end

function A₁projection(f::ModeAnsatz, G::PointGroup)
	return [ModeAnsatz(O*f.q) for O in G.g] |> sum
end

end
