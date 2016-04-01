module Parameters

export RectangleGeometry

type RectangleGeometry{Dim}
	N::NTuple{Dim, Integer}
	Δx⃗::NTuple{Dim, Real} 
	L⃗::NTuple{Dim, Real}
end

function RectangleGeometry{Dim}(N::NTuple{Dim, Real}, Δx⃗::NTuple{Dim, Real})
	L⃗ = ntuple( i -> N[i]*Δx⃗[i], Dim)
	RectangleGeometry(N, Δx⃗, L⃗)
end

function RectangleGeometry{Dim}(dim::Dim, Nₓ::Integer, Δx::Real)
	N = ntuple( i -> Nₓ, Dim)
	Δx⃗= ntuple( i -> Δx, Dim)
	L⃗ = ntuple( i -> Nₓ*Δx, Dim)			
	RectangleGeometry{Dim}(N, Δx⃗, L⃗)
end



end
