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

module BravaisLattices

export BravaisLattice, reciprocal

type BravaisLattice{N}
	a::NTuple{N, Vector}
	function BravaisLattice(a...)
		@assert validbasis(a)
		if N == length(a)
		else error("Wrong basis dimensions for bravais lattice")
		end
		new(a)
	end
end

function reciprocal(B::BravaisLattice{2})
	# Construct reciprocal lattic
	R = [[0 -1]
		 [1  0]]	# 90 degree rotation operator
	b1 = 2π*(-R*B.a[2])/(B.a[1]⋅(-R*B.a[2]))
	b2 = 2π*(-R*B.a[1])/(B.a[2]⋅(-R*B.a[1]))
	return BravaisLattice{2}(b1, b2)
end

function reciprocal(B::BravaisLattice{3})
	b1 = 2π*(B.a[2]×B.a[3])/(B.a[1]⋅(B.a[2]×B.a[3]))
	b2 = 2π*(B.a[3]×B.a[1])/(B.a[2]⋅(B.a[3]×B.a[1]))
	b3 = 2π*(B.a[1]×B.a[2])/(B.a[3]⋅(B.a[1]×B.a[2]))
	return BravaisLattice{3}(b1, b2, b3)
end

function validbasis(a)
	if eltype(a) <: Vector
	else error("All basis vectors must be real vectors of the same length")
	end
	if [length(i) == length(a) for i in a] |> prod
	else error("Inconsistent dimension for basis")
	end
	return true
end

end
