
"""    Stacks with a statically allocated maximum size."""
module StaticStacks
using StaticArrays

mutable struct StaticStack{N,T} <: AbstractVector{T}
	v::MVector{N,T}
	l::Int
	@inline StaticStack{N,T}() where{N,T} = new{N,T}(MVector{N,T}(undef), 0)
end
capacity(::StaticStack{N}) where{N} = N
Base.size(s::StaticStack) = (s.l,)
# Base.checkbounds(::Type{Bool}, s::StaticStack, i) = i ∈ 1:s.l
@inline Base.getindex(s::StaticStack, i) = getindex(s.v, i)
@inline Base.setindex!(s::StaticStack, x, i) = setindex!(s.v, x, i)
@inline function Base.push!(s::StaticStack, x...)
# 	@boundscheck @assert length(s) + length(x) ≤ capacity(s)
	@inbounds s.v[length(s)+1:length(s)+length(x)].= x
	s.l+= length(x)
	return s
end
@inline function Base.pop!(s::StaticStack)
# 	@boundscheck @assert !isempty(s)
	s.l-= 1
	return @inbounds s.v[length(s)+1]
end
function Base.resize!(s::StaticStack, n)
# 	@boundscheck n ∈ 1:length(s.v)
	s.l = n
	return s
end

export StaticStack
end # module
module Trees
# using ResumableFunctions
using StaticArrays
using LazySets: convex_hull
using Printf
using ..StaticStacks
using TOML
const A000669 = (1,1, 2, 5, 12, 33, 90, 261, 766, 2312, 7068, 21965,
68954, 218751, 699534, 2253676, 7305788, 23816743, 78023602, 256738751,
848152864, 2811996972, 9353366564, 31204088381, 104384620070,
350064856815, 1176693361956, 3963752002320)

"""    A number fixed at compile time."""
struct StaticInt{N} <: Signed end
@inline StaticInt(N::Int) = StaticInt{N}()
Base.show(io::IO, ::StaticInt{N}) where{N} = print(io, "StaticInt(", N, ")")
Base.promote_rule(T::Type{<:Signed}, ::Type{<:StaticInt}) = T
Base.convert(T::Type{<:Signed}, ::StaticInt{N}) where{N} = convert(T, N)
Base.convert(::Type{T}, ::T) where{T<:StaticInt} = T()
Base.convert(::Type{<:StaticInt}, ::StaticInt) = error("incompatible values")
Base.iszero(::StaticInt{N}) where{N} = iszero(N)

"""    AbstractTree - a base type describing the interface for costs
attached to trees:
 - `parent_cost(t)`
 - `sibling_cost(t)`
"""
abstract type AbstractTree end

""" A general (non-planar) tree, encoded as a list of node valences,
as seen in post-order traversal (Łukasiewicz word).
Type parameters:
 - `T`: vector-like type holding the degrees of the tree
		(either `Vector{Int}` or `SVector{2*N-1,Int}`);
 - `P`: vector-like type used for paths in the tree
    (either `Vector{Int}` or `StaticStack{N,Int}`).
"""
struct PostTree{V} <: AbstractTree
	degrees::V
end

function PostTree(s::AbstractString)
	v = Int[]
	stack = Int[]
	for c in reverse(s)
		# if c is not a parenthesis then it is treated as a pair "()" of
		# parentheses; this allows dots or hyphens as shorthand
		if c != '(' # opening parenthesis
			push!(stack, 0)
		end
		if c != ')'
			n = pop!(stack)
			for _ in 3:n; push!(v, 1) end
			push!(v, n)
			!isempty(stack) && (stack[end]+= 1)
		end
	end
	return PostTree(v)
end

empty_stack(::Vector) = Int[]
empty_stack(::StaticVector{N}) where{N} = StaticStack{N,Int}()

# PostTree(deg::AbstractVector) =
# 	compute_costs!(PostTree(deg, similar(deg, Int), similar(deg, Int)))
"""    Computes all costs attached to this tree."""
function costs(t::PostTree)
	nl, pc, sc = (empty_stack(t.degrees) for _ in 1:3) # stacks for computing costs
	parent_cost = sibling_cost = n = 0
	for (i, d) in pairs(t.degrees)
		if iszero(d)
			n = 1
			parent_cost = 0
			sibling_cost = 0
		elseif isone(d)
			n = last(nl)
			parent_cost = last(pc)
			sibling_cost = last(sc)
		else
			n = sum(@view nl[end-d+1:end])
			parent_cost = sum(@view pc[end-d+1:end]) + n
			sibling_cost = sum(@view sc[end-d+1:end]) + n*(d-1)
		end
# 		@debug("t[$i]=$d: pc=$pc, new=$parent_cost")
		s = length(nl) - d + 1
		resize!(nl, s); nl[end] = n
# 		@debug("  => now nl=$nl")
		resize!(pc, s); pc[end] = parent_cost
		resize!(sc, s); sc[end] = sibling_cost
	end
	parent_cost, sibling_cost
end
nleaves(t::PostTree) = count(iszero, t.degrees)

@inline Base.copy(t::PostTree) = deepcopy(t)

function Base.show(io::IO, t::PostTree)
	s = empty_stack(t.degrees)
	for d in Iterators.reverse(t.degrees)
		isone(d) && continue # ignore the “filler” nodes
		if iszero(d)
			printstyled(io, '.', color=length(s)%8)
			k = length(s)
			while k in eachindex(s) && isone(s[k])
				printstyled(io, ')', color=(k-1)%8)
				k-=1
			end
			iszero(k) && return
			resize!(s, k); s[k]-= 1
		else
			printstyled(io, '(', color=length(s)%8)
			push!(s, d)
		end
	end
end
function Base.display(t::PostTree)
	println(t)
	level = similar(t.degrees); h = 0
	s = empty_stack(t.degrees)
	for (j,d) in Iterators.reverse(pairs(t.degrees))
		level[j] = length(s)
		isone(d) && continue
		if iszero(d)
			k = findlast(!isone, s)
			isnothing(k) && break
			resize!(s, k); s[k]-= 1
		else
			push!(s, d)
			length(s) > h && (h = length(s))
		end
	end
	for i in 1:length(t.degrees)
		@printf("%3d│%s %d\n", i, "  "^(h-level[i]), t.degrees[i])
	end
end
"""    Writes in the stack a branch in the shape of the first tree of
weight `a`, i.e (((..).).) etc; the post-form is 0(02)^(a-1);
at the `(2a-1)` indices `i:i+2a-2`. Returns the index just after this
tree.  """
@inline function write_firsttree!(t::PostTree, start, a)
	iszero(a) && return start
	t.degrees[start] = 0
	# FIXME: should we also write the intermediate costs?
	for n in 1:a-1
		t.degrees[start+2*n-1] = 0
		t.degrees[start+2*n] = 2
	end
# 	t.parent_cost[start+2*a-2] = t.sibling_cost[start+2*a-2] = (a^2+a-2)÷2
	return start+2*a-1
end
@inline tree_storage(n::Int) = Vector{Int}(undef, 2*n-1)
@inline tree_storage(::StaticInt{N}) where{N} = MVector{2*N-1,Int}(undef)
"""    Returns the first tree with `n` leaves."""
function firsttree(n)
	t = PostTree(tree_storage(n))
	write_firsttree!(t, 1, n)
	return t
end

"""    Returns the first index encoding the branch with given `root`,
i.e. this branch is encoded as indices `start:root`."""
function branch_start(t::PostTree, root)
	s = 0
	for j in root:-1:1
		s+= t.degrees[j]-1
# 		@debug("at t.degrees[$j]=$(t.degrees[j]): s=$s")
		s < 0 && return j
	end
	error("incomplete tree")
end
"""    Returns a pair `(j, v)`, where `j` is the left-most index in this
branch and `v` a vector of all weights of children of this node."""
@inline function children_weights(t::PostTree, parent)
	weights = Int[]; s = 0
	for j in parent:-1:1
		d = t.degrees[j]
		isone(d) && continue
		if s + length(weights) == t.degrees[parent]-1
			push!(weights, 0)
		end
		s+= d-1
		iszero(d) && !isempty(weights) && (weights[end]+= 1)
		s < 0 && return j, reverse!(weights)
	end
end
"""    Finds the (pre-order) last grandparent node,
i.e. the last node of the subtree `t'` formed by
all the “grandparent” nodes, i.e. nodes having at least one
sub-sub-node.

Returns either the path to this grandparent node
(as a stack from root node to grandparent),
or `nothing` if no grandparent node can be found.
"""
function grandparent_node(t::PostTree)
	# `path`: stack of all roots from current node to root of tree
	# `position`: remaining nodes on each element of `path`.
	path, position, nleaves = (empty_stack(t.degrees) for _ in 1:3)
	for (j, d) in Iterators.reverse(pairs(t.degrees))
		isone(d) && continue
		if iszero(d)
			k = something(findlast(!isone, position), 0)
			if k <= length(position) - 2
# 				@debug "pop grandparent: $path at t[$j]=$d; returning $(path[end-1])"
				pop!(path)
				return path
			end
# 			if length(position) - k >= 2
# 				@debug("correct position seems to be $j; path=$path; Returning ($j, $(path[end]), $(path[end-1]))")
# 			end
			if iszero(k) # end of tree
# 				@debug("no grandparent node found in this tree")
				return nothing
			else
				resize!(path, k); resize!(position, k)
				position[k]-= 1
# 				@debug("  "^length(path),"popping, now path, position = $path, $position; $child->$parent")
# 				if !isnothing(parent) && path[k] > child
# 					break
# 				end
			end
		else
# 			child, parent = j, child
			push!(path, j); push!(position, d)
		end
	end
end
""""    Computes next partition of same total,
in rev-lex ordering."""
@inline function nextpartition!(weights)
	n = sum(weights)
	j = findlast(!isone, weights) # there are (L-j) ones to redistribute
	isnothing(j) && return nothing
	w = weights[j]
	# we split (L-j+w) in parts of size (w-1)
	q, r = divrem(length(weights)-j+w, w-1)
	resize!(weights, j+q-iszero(r))
	weights[j:j+q-1].= w-1
	!iszero(r) && (weights[end] = r)
	return weights
end
"""    Advances the subtree from index `root` to the next subtree.
Returns `true` if advancing was successful, `false` otherwise."""
function nexttree!(t::PostTree)
	path = grandparent_node(t)
# 	println("[7m$(t.degrees)[m")
# 	println("grandparent = $path")
	isnothing(path) && return
	parent = last(path)
	### first modify the branch starting at this parent node:
	# compute the list of children weights of this node:
	start, weights = children_weights(t, parent)
# 	println("weights $weights ")
	weights = nextpartition!(weights)
# println("   => $weights")
	isnothing(weights) && return nothing
	self_low = start
	for w in weights
		start = write_firsttree!(t, start, w)
	end
	# fill with ones to avoid moving tree branches around
	t.degrees[start:parent-1].= 1
	t.degrees[parent] = length(weights)
	### then modify all the (rightward) uncles:
	# - those with the same weight as the respective parent should be made
	#   similar to the corresponding parent tree,
	# - those with a different weight (always lower) should be
	#   reinitialized to the first tree (caterpillar) with this weight
	# TODO: `grandparent_node` should return a list of those uncles
	# (index and level in the tree).

	# visit and reset uncles to the right of `parent`,
	# always while preserving node weight:
	self = parent
	self_weight = sum(weights)
	for parent in Iterators.reverse(@view path[1:end-1])
		h = 0
		uncle = 0; nleaves = 0
		n_uncles = 0; uncles_weight = 0
# 		println("  examining interval ]$self:$parent[; weight of $self is $self_weight")
		for j in parent-1:-1:self+1
			d = t.degrees[j]
			isone(d) && continue
# 			println("    at ($j,$d): h=$h")
			if iszero(h)
				uncle = j; nleaves = 0
				h+= 1
			end
			h+= d-1; nleaves+= iszero(d)
			if iszero(h)
# 				println("reset uncle at $j:$uncle to $nleaves")
				if nleaves == self_weight
# 					println("!!! FOUND AN UNCLE WITH THE SAME WEIGHT ($nleaves) AS SELF")
# 					println("uncle is at $j:$uncle, we are at $self_low:$self")
# 					println("$(t.degrees[j:uncle])  ... $(t.degrees[self_low:self])")
					@assert uncle-j == self-self_low
					t.degrees[j:uncle].= @view t.degrees[self_low:self]
				else
					write_firsttree!(t, j, nleaves) # writes j to j+2n-2
					t.degrees[j+2*nleaves-1:uncle].= 1
				end
				n_uncles+= 1
				uncles_weight+= nleaves
			end
# 			h < 0 && break
		end
		# we need to keep track of the current weight of `self`,
		# which may contain leaves added to the *left*.
		# We know that `parent` has `t.degrees[parent]` children,
		# including `self` and the already visited right-siblings.
		# We now explore the left-siblings:
		rem_siblings = t.degrees[parent]-n_uncles-1
# 		println("enlarging left side of current tree:")
		while rem_siblings > 0
			c = 1
			while self_low > 0
				d = t.degrees[self_low]
				iszero(d) && (self_weight+= 1)
				c+= d-1
				if c < 0
					break
				end
				self_low-=1
			end
			rem_siblings-= 1
# 			println("found a sibling at $self_low, $rem_siblings remaining")
		end
		self_weight+= uncles_weight
# 		println("  .now self_weight=$self_weight")
		self = parent
	end
	return t
end
@inline nexttree(t) = nexttree!(deepcopy(t))

"""    An in-place iterator over all trees of a given number of leaves.

The number of such trees is OEIS A000669;
The generating function f verifies the recursive definition
2f = z - 1 + Exp(f), where f is the Pólya exponential.

a_n ~ A⋅B^n/n^(3/2) where B ≈ 3.56
"""
struct AllTreesMut
	nleaves::Int
end
function Base.iterate(a::AllTreesMut, state=nothing)
	if isnothing(state)
		t = firsttree(a.nleaves)
		return t, t
	end
	r = nexttree!(state)
	isnothing(r) && return
	return state, state
end
struct AllTrees{N<:Signed}
	nleaves::N
end
function Base.iterate(a::AllTrees, state=nothing)
	if isnothing(state)
		t = firsttree(a.nleaves)
		return deepcopy(t), t
	end
	r = nexttree!(state)
	isnothing(r) && return
	return deepcopy(state), state
end
# ntrees(n) = n <= 1 ? 1 : sum(prod(ntrees(q) for q in decomp)
# 		for decomp in partitions(n) if length(decomp) > 1)
# Base.IteratorSize(::Type{<:AllTrees}) = Base.HasLength()
Base.length(t::AllTrees) = A000669[t.nleaves]
Base.IteratorEltype(::Type{<:AllTrees}) = Base.HasEltype()
Base.eltype(::AllTrees{T}) where{T<:Signed} = 
	PostTree{Base.return_types(tree_storage,(T,))[1]}

"""    The parent and sibling cost of a tree,
abstracted as a 2-dimensional vector (for taking convex hull)."""
struct TreeCost{T<:AbstractTree} <: AbstractVector{Int}
	tree::T
	parent_cost::Int
	sibling_cost::Int
end
@inline TreeCost(T::AbstractTree) = TreeCost(T, costs(T)...)
Base.size(t::TreeCost) = (2,)
Base.getindex(t::TreeCost, i) = (i == 1) ? t.parent_cost : t.sibling_cost
Base.IndexStyle(::Type{TreeCost}) = IndexLinear()

"""    Computes the best trees with `n` leaves.
Optional parameter `bound` indicates the frequency of convex hull
computation: convex hull is taken whenever the current set
contains more than `bound` trees.
Default value is 512, which empiricaly seems quite good now."""
function compute_best_trees(n; bound=512)
	# first include three dummy points to compute the south-west envelope:
	# note: `n^2` is used here as a majorant of both the pc and the sc of
	# trees; the value nleaves==-1 is a marker for later identifying this
	# particular dummy point (and removing it).
	envelope = TreeCost[]
	x1, x2 = y1, y2 = typemax(Int), typemin(Int)
	c = 0
	for t in AllTrees(n)
		c+=1
		if iszero(c%10^4)
			u = A000669[n]
			@printf("n=%d tree %d/%d (%.2f %%)\r", n, c, u, 100*c/u)
		end
		p = TreeCost(t)
		x1 = min(x1, p[1]); x2 = max(x2, p[1])
		y1 = min(y1, p[2]); y2 = max(y2, p[2])
		push!(envelope, p)
		l = length(envelope)
		if l >= bound && ispow2(l)
# 			println("(weight ", n, ")  with ", l, " points, taking convex hull..")
			envelope = convex_hull(envelope)
# 			println("  => reduced to ", length(envelope))
		end
		for (i, q) in pairs(envelope)
			if (@view p[1:2]) == (@view q[1:2])
				insert!(envelope, i, p)
			end
		end
	end
	push!(envelope,
		TreeCost(PostTree(Int[]), x1, y2+1),
		TreeCost(PostTree(Int[]), x2+1, y2+1),
		TreeCost(PostTree(Int[]), x2+1, y1))
	envelope = convex_hull(envelope)
	circshift!(envelope, -1-findfirst(t->sum(t)==x2+y2+2, envelope))
	resize!(envelope, length(envelope)-3)
	envelope
end
const BEST_TREES=Dict{Int,Vector{TreeCost}}()
const BEST_TREES_FILE="best_trees.toml"

"""    Returns the list of best trees with `n` leaves,
either using the cached value in `BEST_TREES` or computing it anew."""
function best_trees(n; kw...)
	haskey(BEST_TREES, n) && return BEST_TREES[n]
	BEST_TREES[n] = compute_best_trees(n; kw...)
	write_best_trees(BEST_TREES_FILE)
	return BEST_TREES[n]
end
@inline best_trees(range::AbstractUnitRange; kw...) =
	[best_trees(x; kw...) for x ∈ range if x ≥ 3]

"""    Saves all currently known best trees to disk cache."""
function write_best_trees(file)
	open(file, "w") do io
		TOML.print(io, Dict(string(k) =>
			[ Dict("pc" => x.parent_cost, "sc" => x.sibling_cost,
				"tree" => x.tree.degrees) for x in v]
				for (k,v) in pairs(BEST_TREES)))
	end
end
"""    Reads all the previously computed best trees from the cache."""
function read_best_trees(file)
	for (k,v) in pairs(TOML.parsefile(file))
		n = parse(Int, k)
		n < 3 && continue
		BEST_TREES[n] =
			[ TreeCost(PostTree(x["tree"]), x["pc"], x["sc"]) for x in v ]
	end
end

function __init__()
	run(`sh -c 'rm -f /tmp/jl_\*'`)
	isfile(BEST_TREES_FILE) && read_best_trees(BEST_TREES_FILE)
end

"""    Outputs a single `PostTree` to an IO as a LaTeX document,
using TikZ and TikZ-qtree packages."""
function tikz(io::IO, t::PostTree, pc=0, sc=0)
	print(io, "\$(", pc, ", ", sc, ")\$:\n\n")
	print(io, "\\Tree")
	s = empty_stack(t.degrees)
	for d in Iterators.reverse(t.degrees)
		isone(d) && continue
		if iszero(d)
			print(io, "{} ")
			k = length(s)
			while k in eachindex(s) && isone(s[k])
				print(io, "] ")
				k-= 1
			end
			iszero(k) && break
			resize!(s, k); s[k]-= 1
		else
			print(io, "[.{} ")
			push!(s, d)
		end
	end
end
function tikz(io::IO, v::AbstractVector{<:TreeCost})
	tikz(io, v[1].tree, v[1][1], v[1][2])
	for i in 2:length(v)
		((x1, y1), (x2, y2)) = v[i-1:i]
		λ = (y2-y1)/(x2-x1)
		println(io, "\n\n(\$\\lambda = ", λ, "\$)\n\n")
		println(io, "\\hbox to \\hsize{%")
		tikz(io, v[i].tree, v[i][1], v[i][2])
		println(io, "\\hfill")
		f=tempname()
		graphviz(v[i].tree, f)
		println(io, "\n\\vtop{\\vss\\includegraphics[scale=.30]{$f.pdf}}\\hfill}\n\n")
	end
end
""""    Merges all nodes of the same degree, returning the reduced graph
(as a dictionary containing the adjacency matrix)."""
function reduced_graph(t::PostTree)
	links = Dict{Int,Set{Int}}()
	path = empty_stack(t.degrees) # path from root to current node
	s = empty_stack(t.degrees) # remaining siblings of current node
# 	push!(path, -1)
	for d in Iterators.reverse(t.degrees)
		isone(d) && continue
		!isempty(path) && push!(get!(links, last(path), Set{Int}()), d)
		if iszero(d)
			k = findlast(!isone, s)
			isnothing(k) && break
			resize!(s, k); resize!(path, k); s[k]-= 1
		else
			push!(s, d)
			push!(path, d)
		end
	end
	return links
end
"""    Outputs a reduced graph as a Graphviz `.dot` file."""
function graphviz(links::Dict, filename=tempname())
	open(filename*".gv", "w") do io
		println(io, """digraph G {
	nodesep=0.1;
	pad=0;
	node[fontsize=24,shape=box, margin=0, height=0, width=0];
	bgcolor="transparent";
""")
		node = Dict{Int,String}()
# 		for d in union(values(links)..., keys(links))
# 			if d < 0
# 				node[d] = "init"
# 				println(io, node[d], "[label=\"\", shape=plaintext];")
# 			elseif iszero(d)
# 				node[d] = "leaf"
# 				println(io, node[d], "[label=\"\", shape=plaintext];")
# 			else
# 				node[d] = "node$d"
# 				println(io, node[d], "[label=\"", d, "\", shape=box, style=filled, color=lightgrey];")
# 			end
# 		end
		for d1 in keys(links)
			c = (0 ∈ links[d1]) ? "00cc88" : "f8f8f8"
			println(io, "node$d1[label=\"$d1\", color=\"#$c\"];")
			for d2 in links[d1]
				iszero(d2) && continue
				println(io, "node$d1 -> node$d2;")
			end
		end
		println(io, "}")
	end
	cd(dirname(filename)) do
		run(`dot -Tpdf $filename.gv -o $filename.pdf`)
	end
	return filename
end
@inline graphviz(tree::PostTree, fn...) = graphviz(reduced_graph(tree), fn...)
@inline graphviz(tree::TreeCost, fn...) = graphviz(tree.tree, fn...)
function report(file="/tmp/trees_report.tex")
	open(file, "w") do io
		println(io, raw"""
\documentclass{article}
\usepackage{tikz}
\usepackage{tikz-qtree}
\usepackage{datetime}
\usepackage{graphicx}
\begin{document}
Compiled \today\ \currenttime """)
		for (k,v) in sort(pairs(BEST_TREES); by=first)
			println("\n\e[7m Trees with ", k, " leaves \e[m")
			println(io, "\n\\newpage{\\Large\\textbf{Best ($(length(v))) trees with \$$k\$ leaves}}\n")
			tikz(io, v)
		end
		println(io, "\\end{document}")
	end
	cd(dirname(file)) do
		run(`pdflatex $(basename(file))`)
	end
end
function graph(file="/tmp/trees_graph.tex")
	open(file, "w") do io
		println(io, raw"""
\documentclass{article}
\usepackage{tikz}
\usepackage{graphicx}
\usepackage{datetime}
\usepackage[margin=15mm,paperwidth=40cm,paperheight=40cm]{geometry}
\pagestyle{empty}
\begin{document}
Compiled \today\ \currenttime

\begin{tikzpicture}[scale=1.8]""")
		for (n,v) in sort(pairs(BEST_TREES); by=first)
			λ = Float64[]
			y = n+.5
			for i in 2:length(v)
				((x1, y1), (x2, y2)) = v[i-1:i]
				push!(λ, -(y2-y1)/(x2-x1))
			end
			println(io, """

% n=$n, λ=$λ
\\node at (-.2,$y) {\$$n\$};
\\draw[very thick] (0,$n)--(0,$(n+1));
\\draw (0,$n)--($(maximum(λ)+1),$n);
\\draw (0,$(n+1))--($(maximum(λ)+1),$(n+1));
""")
			for m in λ
			println(io, "\\draw[thin] ($m, $n)--($m, $(n+1));")
			end
			for (i, t) in pairs(v)
				f = graphviz(t.tree)
				x = (get(λ, i-1, maximum(λ)+1) + get(λ, i, 0))/2
				y1 = y - .1 + .2*(i%2)
				println(io, """
% (x=$x) best tree $i: $(t.tree.degrees)
\\node at ($x, $y1) {\\includegraphics[scale=.30]{$f.pdf}};""")
			end
		end
		println(io, raw"""
\end{tikzpicture}
\end{document}""")
	end
	run(`cat $file`)
	cd(dirname(file)) do
		run(`pdflatex $(basename(file))`)
	end
end
end # module

function test_vectors(filename::AbstractString)
	open(filename) do f
		for line in eachline(f)
			values = split(line, r"\s*/\s*")
			length(values) == 6 || continue
			nleaves = tryparse(Int, values[1])
			isnothing(nleaves) && continue
			pc, sc = parse.(Int, values[4:5])
			tree = Trees.PostTree(values[6])
			v6 = replace(values[6], '-' => '.')
			if v6 ≠ repr(tree)
				println("[41mBAD TREE[m: $(values[6])\n$v6\n$(tree.degrees)\n$tree")
				return
			end
			pc1, sc1 = Trees.costs(tree)
			if (pc1, sc1) ≠ (pc, sc)
				println("[43mBAD COSTS[m: $pc1, $sc1, tabulated $pc, $sc")
				return
			end
		end
	end
end
# t = Trees.PostTree("((.(..))(..))")
# 
# t0=Trees.PostTree([0,0,2,0,2,0,0,0,1,3,2])
# t0=Trees.PostTree([0, 0, 2, 0, 2, 0, 2, 0, 0, 0, 1, 3, 2, 0, 2, 0, 2])
# t0=Trees.PostTree([0,0,2,0,2,0,0,0,1,3,2])
# t1=Trees.nexttree(t0)
# t1=Trees.nexttree!(deepcopy(t0))
t=Trees.PostTree([0, 0, 0, 1, 3, 0, 0, 0, 1, 3, 0, 0, 1, 1, 4])
run(`sh -c 'rm -f /tmp/jl_\* /tmp/nodes*'`)
# Trees.graphviz(t, "/tmp/nodes")
