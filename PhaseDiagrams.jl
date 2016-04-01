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

using Functional, PointGroups, BravaisLattices
export ModeAnsatz, A₁projection, mode


"""
	ModeAnsatz

Make single mode ansatz from a reciprocal lattice vector q.

Given a reciprocal lattice vector q, the ansatz function is of the form:

	\$\$ f(x) = cos(q⋅x) \$\$

The single mode ansatz is callable, as it is a function. The call
signature for a mode ansatz f is
	f(x::Vector)

**Summary:**
	ModeAnsatz <: Function
**Subtypes:**
	None
"""
type ModeAnsatz <: Function
	q::Vector
	F::Function
	function ModeAnsatz(q::Vector)
		f = (x::Vector) -> cos(q⋅x)
		new(q, f)
	end
end

"""
	ModeAnsatz(q::BravaisLattice)

Create mode ansatz from the basis of a bravais lattice.

The bravais lattice should be the reciprocal lattice of some lattice of
interest. If this mode ansatz is passed to the A₁projection the result will
be a one-mode approximation.
"""
function ModeAnsatz(q::BravaisLattice)
	[ModeAnsatz(k) for k in q.basis] |> sum
end

"""
	mode(f::ModeAnsatz)

Return reciprocal lattice mode of ansatz
"""
function mode(f::ModeAnsatz)
	return f.q
end

function (f::ModeAnsatz)(x::Vector)
	return f.F(x)
end

"""
	A₁projection(f::ModeAnsatz, G::PointGroup)

Return A₁ symmetry adapted ansatz of f for the group G
"""
function A₁projection(f::ModeAnsatz, G::PointGroup)
	return [ModeAnsatz(O*f.q) for O in G.g] |> sum
end

end
