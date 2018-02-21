using DataStructures

function filter_components(input, min_size)
    if length(size(input)) == 2
        return filter_components2d(input, min_size)
    elseif length(size(input)) == 3
        return filter_components3d(input, min_size)
    end
end

function filter_components3d(input, min_size)
    output = zeros(Int64, size(input))
    (n, m, d) = size(input)
    sets = IntDisjointSets(0)
    sizes = counter(Int64)
    for k in 1:d, j in 1:m, i in 1:n
        if !is_fg(input[i, j, k])
            continue
        end
        components = neighbor_components3d(i, j, k, output)
        if length(components) == 0
            component = push!(sets)
            output[i, j, k] = component
            push!(sizes, component)
        else
            output[i, j, k] = components[1]
            push!(sizes, components[1])
            if length(components) == 2
                union!(sets, components...)
            end
        end
    end

    # avoid modifying counter while iterating over it
    size_list = [(key, count) for (key, count) in sizes]
    for (key, count) in size_list
        root = find_root(sets, key)
        push!(sizes, root, reset!(sizes, key))
    end
    for k in 1:d, j in 1:m, i in 1:n
        if !is_fg(input[i, j, k])
            continue
        end
        component = output[i, j, k]
        component_size = sizes[find_root(sets, component)]
        if component_size < min_size
            output[i, j, k] = 0
        else
            output[i, j, k] = 1
        end
    end
    return output
end

function filter_components2d(input, min_size)
    output = zeros(Int64, size(input))
    (n, m) = size(input)
    sets = IntDisjointSets(0)
    sizes = counter(Int64)
    for j in 1:m, i in 1:n
        if !is_fg(input[i, j])
            continue
        end
        components = neighbor_components2d(i, j, output)
        if length(components) == 0
            component = push!(sets)
            output[i, j] = component
            push!(sizes, component)
        else
            output[i, j] = components[1]
            push!(sizes, components[1])
            if length(components) == 2
                union!(sets, components...)
            end
        end
    end

    # avoid modifying counter while iterating over it
    size_list = [(key, count) for (key, count) in sizes]
    for (key, count) in size_list
        root = find_root(sets, key)
        push!(sizes, root, reset!(sizes, key))
    end
    for j in 1:m, i in 1:n
        if !is_fg(input[i, j])
            continue
        end
        component = output[i, j]
        component_size = sizes[find_root(sets, component)]
        if component_size < min_size
            output[i, j] = 0
        else
            # output[i, j] = find_root(sets, component)
            output[i, j] = 1
        end
    end
    return output
end

is_fg(x) = x == one(x)

function neighbor_components2d(i, j, labels)
    if i > 1 && j > 1
        neighbors = [(i-1, j), (i, j-1)]
    elseif i > 1
        neighbors = [(i-1, j)]
    elseif j > 1
        neighbors = [(i, j-1)]
    else
        neighbors = Array{Int64}(0)
    end
    components = filter(x -> x > 0, [labels[idx...] for idx in neighbors])
    sort!(components)
end

function neighbor_components3d(i, j, k, labels)
    if i > 1 && j > 1 && k > 1
        neighbors = [(i-1, j, k), (i, j-1, k), (i, j, k-1)]
    elseif i > 1 && j > 1
        neighbors = [(i-1, j, k), (i, j-1, k)]
    elseif i > 1 && k > 1
        neighbors = [(i-1, j, k), (i, j, k-1)]
    elseif j > 1 && k > 1
        neighbors = [(i, j-1, k), (i, j, k-1)]
    elseif i > 1
        neighbors = [(i-1, j, k)]
    elseif j > 1
        neighbors = [(i, j-1, k)]
    elseif k > 1
        neighbors = [(i, j, k-1)]
    else
        neighbors = Array{Int64}(0)
    end
    components = filter(x -> x > 0, [labels[idx...] for idx in neighbors])
    sort!(components)
end
