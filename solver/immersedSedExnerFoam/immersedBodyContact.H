/*---------------------------------------------------------------------------*\
  Parallel-safe multi-point Winkler/Coulomb contact with the Exner bed.
\*---------------------------------------------------------------------------*/

contactForce = vector
(
    0,
    0,
    contactFluidDensity*displacedVolume*mag(g.value())
);
contactMoment = Zero;

const vector bedNormal(0, 0, 1);
const vectorField& bedCentres = bed.aMesh().areaCentres();

// Project only downward-facing surface patches onto the horizontal bed.  The
// summed projected area is 0.01 m2 upright and 0.02 m2 on either long side.
scalarField currentPointAreas(refContactPoints.size(), Zero);
scalar currentProjectedArea = 0;

forAll(refContactPoints, pointI)
{
    if (useDirectionalContactQuadrature)
    {
        const vector currentNormal =
            bodyMotion.orientation() & refContactNormals[pointI];
        currentPointAreas[pointI] =
            refContactPointAreas[pointI]
           *max(-(currentNormal & bedNormal), scalar(0));
    }
    else
    {
        currentPointAreas[pointI] =
            contactArea/scalar(refContactPoints.size());
    }
    currentProjectedArea += currentPointAreas[pointI];
}

forAll(refContactPoints, pointI)
{
    const scalar pointArea = currentPointAreas[pointI];
    if (pointArea <= SMALL)
    {
        continue;
    }

    const scalar pointMass =
        bodyMotion.mass()*pointArea/(currentProjectedArea + VSMALL);
    const scalar pointStiffness = foundationStiffness*pointArea;
    const scalar pointDamping =
        2.0*contactDampingRatio
       *Foam::sqrt(doubleScalar(max(pointStiffness*pointMass, SMALL)));
    const scalar maxPointForce = bearingCapacity*pointArea;

    const point contactPoint = bodyMotion.transform(refContactPoints[pointI]);
    const vector contactVelocity = bodyMotion.velocity(contactPoint);

    scalar localDistanceSqr = GREAT;
    scalar localBedZ = 0;

    forAll(bedCentres, bedFaceI)
    {
        const scalar dx = contactPoint.x() - bedCentres[bedFaceI].x();
        const scalar dy = contactPoint.y() - bedCentres[bedFaceI].y();
        const scalar distanceSqr = sqr(dx) + sqr(dy);

        if (distanceSqr < localDistanceSqr)
        {
            localDistanceSqr = distanceSqr;
            localBedZ = bedCentres[bedFaceI].z();
        }
    }

    scalar globalDistanceSqr = localDistanceSqr;
    reduce(globalDistanceSqr, minOp<scalar>());

    const bool ownsNearest =
        mag(localDistanceSqr - globalDistanceSqr)
      < max(SMALL, 1.0e-10*max(globalDistanceSqr, scalar(1)));

    scalar bedZContribution = ownsNearest ? localBedZ : 0;
    label nearestCount = ownsNearest ? 1 : 0;
    reduce(bedZContribution, sumOp<scalar>());
    reduce(nearestCount, sumOp<label>());

    if
    (
        nearestCount == 0
     || globalDistanceSqr > sqr(contactSearchRadius)
    )
    {
        continue;
    }

    const scalar bedZ = bedZContribution/scalar(nearestCount);
    const scalar penetration = bedZ - contactPoint.z();

    if (penetration <= 0)
    {
        continue;
    }

    const scalar normalVelocity = contactVelocity & bedNormal;
    scalar normalForce =
        pointStiffness*penetration - pointDamping*normalVelocity;
    normalForce = min(max(normalForce, scalar(0)), maxPointForce);

    const vector tangentialVelocity =
        contactVelocity - normalVelocity*bedNormal;
    const scalar tangentialSpeed = mag(tangentialVelocity);
    vector tangentialForce(Zero);

    if (tangentialSpeed > SMALL)
    {
        const scalar mu =
            dynamicFriction
          + (staticFriction - dynamicFriction)
           *Foam::exp(-sqr(tangentialSpeed/(frictionVelocity + VSMALL)));

        tangentialForce =
            -mu*normalForce
            *Foam::tanh
            (
                doubleScalar(tangentialSpeed/(frictionVelocity + VSMALL))
            )
            *tangentialVelocity/(tangentialSpeed + VSMALL);
    }

    const vector pointForce = normalForce*bedNormal + tangentialForce;
    contactForce += pointForce;
    contactMoment +=
        (contactPoint - bodyMotion.centreOfMass()) ^ pointForce;
}
