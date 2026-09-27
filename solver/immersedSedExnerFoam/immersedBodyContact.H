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

const scalar pointArea = contactArea/scalar(refContactPoints.size());
const scalar pointMass = bodyMotion.mass()/scalar(refContactPoints.size());
const scalar pointStiffness = foundationStiffness*pointArea;
const scalar pointDamping =
    2.0*contactDampingRatio
   *Foam::sqrt(doubleScalar(max(pointStiffness*pointMass, SMALL)));
const scalar maxPointForce = bearingCapacity*pointArea;
const vector bedNormal(0, 0, 1);
const vectorField& bedCentres = bed.aMesh().areaCentres();

forAll(refContactPoints, pointI)
{
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
