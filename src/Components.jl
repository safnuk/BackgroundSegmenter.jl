using DataStructures

function filter_components(input, min_size)
    output = zeros(Int64, size(input))
    filter_components!(output, input, min_size)
end

function filter_components!(input, min_size)
    filter_components!(input, input, min_size)
end

function filter_components!(output, input, min_size)
    dims = length(size(input))
    sets = IntDisjointSets(0)
    sizes = counter(Int64)
    min_bound = CartesianIndex(ones(Int, dims)...)
    myEye = eye(Int64, dims)
    offsets = [CartesianIndex(myEye[x, :]...) for x in 1:dims]
    for idx in CartesianRange(size(input))
        if !is_fg(input[idx])
            continue
        end
        components = neighbor_components(idx, output, offsets, min_bound)
        if length(components) == 0
            component = push!(sets)
            output[idx] = component
            push!(sizes, component)
        else
            output[idx] = components[1]
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
    labels = Dict{Int, Int}(-1 => 0)
    for idx in CartesianRange(size(input))
        if !is_fg(input[idx])
            continue
        end
        component = output[idx]
        root = find_root(sets, component)
        component_size = sizes[root]
        if component_size < min_size
            output[idx] = 0
        else
            label = get_label!(labels, root)
            output[idx] = label
        end
    end
    return output
end

is_fg(x) = x > zero(x)

function get_label!(labels, key)
    if haskey(labels, key)
        return labels[key]
    else
        labels[-1] += 1
        label = labels[-1]
        labels[key] = label
        return label
    end
end

function neighbor_components(idx, labels, offsets, min_bound)
    components = []
    for offset in offsets
        neighbor = idx - offset
        if max(neighbor, min_bound) == neighbor && labels[neighbor] > 0
            components = vcat(components, [labels[neighbor]])
        end
    end
    components
end
