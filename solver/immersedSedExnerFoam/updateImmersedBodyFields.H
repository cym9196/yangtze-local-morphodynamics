/*---------------------------------------------------------------------------*\
  Update smooth cuboid mask and local rigid-body velocity on the fixed mesh.
\*---------------------------------------------------------------------------*/

{
    const vector bodyCentre = bodyMotion.centreOfMass();
    const tensor bodyQT = bodyMotion.orientation().T();
    const vectorField& cellCentres = mesh.C();

    forAll(cellCentres, cellI)
    {
        const vector localPoint = bodyQT & (cellCentres[cellI] - bodyCentre);
        const vector q
        (
            mag(localPoint.x()) - bodyHalfSize.x(),
            mag(localPoint.y()) - bodyHalfSize.y(),
            mag(localPoint.z()) - bodyHalfSize.z()
        );

        const vector qPositive
        (
            max(q.x(), scalar(0)),
            max(q.y(), scalar(0)),
            max(q.z(), scalar(0))
        );

        const scalar signedDistance =
            mag(qPositive)
          + min(max(q.x(), max(q.y(), q.z())), scalar(0));

        bodyMask[cellI] =
            0.5
           *(
                1.0
              - Foam::tanh
                (
                    doubleScalar
                    (
                        signedDistance/(interfaceThickness + VSMALL)
                    )
                )
            );

        bodyVelocity[cellI] = bodyMotion.velocity(cellCentres[cellI]);
    }

    bodyMask.correctBoundaryConditions();
    bodyVelocity.correctBoundaryConditions();
}
