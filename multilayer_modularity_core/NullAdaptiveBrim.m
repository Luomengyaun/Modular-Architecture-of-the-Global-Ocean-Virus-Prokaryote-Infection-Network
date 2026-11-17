classdef NullAdaptiveBrim < AdaptiveBrim
    % AdaptiveBrim for null statistics; omit only final display sorting.
    % Detect computes Qb and Qr before SortModules. All optimization,
    % component handling and random draws remain in the original class.
    methods
        function obj = NullAdaptiveBrim(matrix)
            obj = obj@AdaptiveBrim(matrix);
        end
    end
    methods(Access = protected)
        function obj = SortModules(obj) %#ok<MANU>
            % No plots or exported null partitions need display ordering.
        end
    end
end
