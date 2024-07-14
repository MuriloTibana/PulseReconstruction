% Return the value of exp(a) * erf(z) using various approximations
function erf = error_function(z,a)

%===== Calculate Using Continued Fraction ===
sub_z = z;
erf = sign(real(sub_z)).*exp(a) - exp(a - sub_z.^2).*sub_z./(sqrt(pi).*(sub_z.^2 + 0.5 ...
        ./(1 + 1./(1.5 + sub_z.^2))));

%===== Patch Poles with Pade Approximation ==
sub_z = z(abs(z+2i) < 1.1);
if not(isempty(sub_z))
    erf(abs(z+2i) < 1.1) = exp(a).*(-18.564802i + 12.326185.*(2i + sub_z) + 13.9957261i.*(2i + sub_z).^2 ...
        - 7.0912286.*(2i + sub_z).^3 - 2.1539997i.*(2i + sub_z).^4 + 0.80057144.*(2i + sub_z).^5)...
        ./ (1 - 2.6545518i.*(2i + sub_z) - 2.9260194.*(2i + sub_z).^2 + 1.665267i.*(2i + sub_z).^3 ...
        + 0.48178231.*(2i + sub_z).^4 - 0.054052386i.*(2i + sub_z).^5);
end

%===== Patch Poles with Pade Approximation ==
sub_z = z(abs(z-2i) < 1.1);
if not(isempty(sub_z))
    erf(abs(z-2i) < 1.1) = exp(a).*(18.564802i + 12.326185.*(-2i + sub_z) - 13.9957261i.*(-2i + sub_z).^2 ...
        - 7.0912286.*(-2i + sub_z).^3 + 2.1539997i.*(-2i + sub_z).^4 + 0.80057144.*(-2i + sub_z).^5)...
        ./ (1 + 2.6545518i.*(-2i + sub_z) - 2.9260194.*(-2i + sub_z).^2 - 1.665267i.*(-2i + sub_z).^3 ...
        + 0.48178231.*(-2i + sub_z).^4 + 0.054052386i.*(-2i + sub_z).^5);
end

%===== Patch Poles with Pade Approximation ==
sub_z = z(abs(z) < 1.4);
if(abs(sub_z) < 1.4)
    erf(abs(z) < 1.4) = sign(real(sub_z)).*sqrt(exp(2.*a) - exp(2.*a - (1.2732395 + 0.14001229.*sub_z.^2) ...
        ./(1 + 0.14001229.*sub_z.^2).*sub_z.^2));
end

end % Function end



