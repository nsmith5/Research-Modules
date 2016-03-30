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
