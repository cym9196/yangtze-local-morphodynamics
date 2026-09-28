/*---------------------------------------------------------------------------*\
  Keep RANS fields inside physically meaningful limits after mesh motion.
\*---------------------------------------------------------------------------*/

volScalarField& boundedK = mesh.lookupObjectRef<volScalarField>("k");
volScalarField& boundedOmega = mesh.lookupObjectRef<volScalarField>("omega");
volScalarField& boundedNut = mesh.lookupObjectRef<volScalarField>("nut");

boundedK.max(dimensionedScalar("kMin", boundedK.dimensions(), kFieldMin));
boundedK.min(dimensionedScalar("kMax", boundedK.dimensions(), kFieldMax));
boundedOmega.max
(
    dimensionedScalar("omegaMin", boundedOmega.dimensions(), omegaFieldMin)
);
boundedOmega.min
(
    dimensionedScalar("omegaMax", boundedOmega.dimensions(), omegaFieldMax)
);
boundedNut.max
(
    dimensionedScalar("nutMin", boundedNut.dimensions(), 0.0)
);
boundedNut.min
(
    dimensionedScalar("nutMax", boundedNut.dimensions(), nutFieldMax)
);

// turbulence->correct() has already evaluated wall functions. Re-running
// their boundary update here would require the temporary kOmegaSST:G field
// after that field has left the registry, so only bound the stored values.
