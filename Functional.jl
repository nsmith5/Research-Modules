module Functional

# Simple module for treating functions algebriacally

import Base: +, -, ^, *, /, sum
export +, -, ^, *, /, sum

*(f::Function, g::Function) 	= (x) -> f(x)*g(x)
+(f::Function, g::Function) 	= (x) -> f(x)+g(x)
-(f::Function) 					= (x) -> -f(x)
-(f::Function, g::Function) 	= (x) -> f(x) - g(x)
^(f::Function, n::Integer) 		= (x) -> f(x)^n
^(f::Function, n::Real)			= (x) -> f(x)^n
^(f::Function, g::Function) 	= (x) -> f(x)^g(x)
/(f::Function, g::Function) 	= (x) -> f(x)/g(x)
/(f::Function, n::Number) 		= (x) -> f(x)/n
/(n::Number, f::Function) 		= (x) -> n/f(x)

sum(F::Array{Function, 1})		= (x) -> sum([A[i](x) for i in F])

end
