module PhaseDiagrams

using Functional
export PointGroup, oblique, rectangular, centered_rectangular, hexagonal, square
export BravaisLattice, reciprocal
export ModeAnsatz, A₁projection, mode

###  Point Groups  ###

type PointGroup
	g::NTuple		# Tuple of matrices
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

### Bravais Lattices ###

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
