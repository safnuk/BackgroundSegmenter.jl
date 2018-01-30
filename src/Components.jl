using DataStructures

function label_components(input, min_size)
    output = zeros(Int64, size(input)...)
    (n, m) = size(input)
    sets = IntDisjointSets(0)
    sizes = counter(Int64)
    for j in 1:m, i in 1:n
        if !is_fg(input[i, j])
            continue
        end
        components = neighbor_components(i, j, output)
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

function neighbor_components(i, j, labels)
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
