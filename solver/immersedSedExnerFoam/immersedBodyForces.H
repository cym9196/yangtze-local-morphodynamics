/*---------------------------------------------------------------------------*\
  Integrate the reaction of the Brinkman forcing and advance the iron block.
\*---------------------------------------------------------------------------*/

hydroForce = Zero;
hydroMoment = Zero;

const vector currentCentre = bodyMotion.centreOfMass();
const scalar invPenaltyTime = 1.0/penaltyTime.value();
const scalar fluidDensityValue = rhoF.value();
const scalarField& cellVolumes = mesh.V();
const vectorField& cellCentresForForce = mesh.C();

forAll(cellVolumes, cellI)
{
    const vector cellForce =
        fluidDensityValue
       *bodyMask[cellI]
       *invPenaltyTime
       *(U[cellI] - bodyVelocity[cellI])
       *cellVolumes[cellI];

    hydroForce += cellForce;
    hydroMoment +=
        (cellCentresForForce[cellI] - currentCentre) ^ cellForce;
}

reduce(hydroForce, sumOp<vector>());
reduce(hydroMoment, sumOp<vector>());

#include "immersedBodyContact.H"

const vector gravityForce = bodyMotion.mass()*g.value();
const vector gravityMoment =
    bodyMotion.mass()*(bodyMotion.momentArm() ^ g.value());

if (runTime.value() >= bodyReleaseTime)
{
    bodyMotion.update
    (
        true,
        hydroForce + gravityForce + contactForce,
        hydroMoment + gravityMoment + contactMoment,
        runTime.deltaTValue(),
        runTime.deltaT0Value()
    );
}

Info<< "immersed body: centre=" << bodyMotion.centreOfMass()
    << " velocity=" << bodyMotion.v()
    << " omega=" << bodyMotion.omega()
    << " force=" << hydroForce
    << " contact=" << contactForce
    << " moment=" << hydroMoment + contactMoment << endl;
