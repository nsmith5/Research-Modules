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

module PointGroups

export PointGroup, oblique, rectangular, centered_rectangular, hexagonal, square

# Group is defined by a Tuple of matrices

type PointGroup
	g::NTuple					# Could be more specific NTuple{N, Matrix...}??
	function PointGroup(g...)
		new(g)
	end
end

# 2D group generators

const E = [[1.0 0.0];
		   [0.0 1.0]]

function Cₙ(n::Integer)
	return [[cospi(2/n) -sinpi(2/n)];
			[sinpi(2/n)  cospi(2/n)]]
end

const σ = [[1  0];
		   [0 -1]]

# All possible 2-D point groups for crystal lattices
const C₂ = PointGroup(E, Cₙ(2))
const C₃ = PointGroup(E, Cₙ(3), Cₙ(3)^2)
const C₄ = PointGroup(E, Cₙ(4), Cₙ(4)^2, Cₙ(4)^3)
const C₆ = PointGroup(E, Cₙ(6), Cₙ(6)^2, Cₙ(6)^3, Cₙ(6)^4, Cₙ(6)^5)
const D₂ = PointGroup(E, Cₙ(2), σ, Cₙ(2)*σ)
const D₃ = PointGroup(E, Cₙ(3), Cₙ(3)^2, σ, Cₙ(3)*σ, Cₙ(3)^2*σ)
const D₄ = PointGroup(E, Cₙ(4), Cₙ(4)^2, Cₙ(4)^3, σ, Cₙ(4)*σ, Cₙ(4)^2*σ, Cₙ(4)^3*σ)
const D₆ = PointGroup(E, Cₙ(6), Cₙ(6)^2, Cₙ(6)^3, Cₙ(6)^4, Cₙ(6)^5, σ,
					   Cₙ(6)*σ, Cₙ(6)^2*σ, Cₙ(6)^3*σ, Cₙ(6)^4*σ, Cₙ(6)^5*σ)

# Alias proper groups to 2D lattice names
const oblique = C₂
const rectangular = D₂
const centered_rectangular = D₂
const hexagonal = D₆
const square = D₄

end
