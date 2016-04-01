module Laplacian

using Parameters
export make_∇²

function index_to_k²(i::Integer, j::Integer, geom::RectangleGeometry{2})
	N = geom.N
	L = geom.L⃗
	kx² = 0.0
	ky² = 0.0
	if i < N[2]>>1
		kx² = (2π*(i-1)/L[1])^2
	else
		kx² = (2π*(i-N-1)/L[1])^2
	end
	if j < N[1]>>1
		ky² = (2π*(j-1)/L[2])^2
	else
		ky² = (2π*(j-N-1)/L[2])^2
	end
	return kx² + ky²	
end

function make_∇²(geom::RectangleGeometry{2})
	N = geom.N
	∇² = Array(Float64, N)
	for j in 1:N[2]
		for i in 1:N[1]
			∇²[i,j] = -index_to_k²(i, j, geom)
		end
	end
	return ∇²
end

end
