function values = apply_phase_operator_local(field_in, phase_type)
%APPLY_PHASE_OPERATOR_LOCAL Local copy of the VWA phase operator helper.

phase_type = lower(char(phase_type));

switch phase_type
    case 'real'
        values = real(field_in);
    case 'imag'
        values = imag(field_in);
    case 'neg_imag'
        values = -imag(field_in);
    otherwise
        error('Unsupported phase operator: %s', phase_type);
end
end
