function Fig_nb = next_fig()
    existing = sort(arrayfun(@(x) x.Number, findobj('type','figure')));
    Fig_nb = existing(find(diff(existing) > 2, 1, 'first')) + 1;  
    if isempty(existing) && isempty(Fig_nb)
        Fig_nb = 1;
    elseif isempty(Fig_nb)
        Fig_nb = existing(end) + 1;
    end
end

