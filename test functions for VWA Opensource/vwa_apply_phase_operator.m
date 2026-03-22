function values = vwa_apply_phase_operator(field_in, phase_type)
%VWA_APPLY_PHASE_OPERATOR  Extract the requested physical phase component.
%
% Supported phase operators:
%   'real'     -> real(field_in)
%   'imag'     -> imag(field_in)
%   'neg_imag' -> -imag(field_in)

    switch lower(phase_type)
        case 'real'
            values = real(field_in);
        case 'imag'
            values = imag(field_in);
        case 'neg_imag'
            values = -imag(field_in);
        otherwise
            error('vwa_apply_phase_operator: unsupported phase operator ''%s''.', phase_type);
    end
end
