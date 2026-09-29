// 中文说明：读取铁块参数，创建平滑浸没掩膜、刚体状态、六自由度运动和接触变量。
/*---------------------------------------------------------------------------*\
  Create the Brinkman immersed cuboid and its six-degree-of-freedom state.
\*---------------------------------------------------------------------------*/

Info<< "Reading immersedBodyProperties\n" << endl;

IOdictionary immersedBodyProperties
(
    IOobject
    (
        "immersedBodyProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const vector bodyHalfSize
(
    immersedBodyProperties.get<vector>("halfSize")
);

const dimensionedScalar penaltyTime
(
    "penaltyTime",
    dimTime,
    immersedBodyProperties.get<scalar>("penaltyTime")
);

const scalar interfaceThickness =
    immersedBodyProperties.getOrDefault<scalar>("interfaceThickness", 0.0025);

const scalar sedimentMaskHalfFraction =
    immersedBodyProperties.getOrDefault<scalar>
    (
        "sedimentMaskHalfFraction",
        0.70
    );

const scalar bodyReleaseTime =
    immersedBodyProperties.getOrDefault<scalar>("releaseTime", 1.0);

// The hydrodynamic field can be spun up while both the rigid body and the
// sediment exchange remain frozen.  Keeping this distinct from releaseTime
// also makes it possible to study a fixed obstacle over a mobile bed.
const scalar sedimentReleaseTime =
    immersedBodyProperties.getOrDefault<scalar>
    (
        "sedimentReleaseTime",
        bodyReleaseTime
    );

const dictionary& contactProperties =
    immersedBodyProperties.subDict("contactProperties");

const pointField refContactPoints
(
    contactProperties.get<pointField>("contactPoints")
);
const scalar contactArea = contactProperties.get<scalar>("contactArea");

// Optional orientation-aware surface quadrature.  Legacy cases that provide
// only contactPoints retain the original equal-area behaviour.
const bool useDirectionalContactQuadrature =
    contactProperties.found("contactNormals")
 && contactProperties.found("contactPointAreas");

vectorField refContactNormals(refContactPoints.size(), vector::zero);
scalarField refContactPointAreas
(
    refContactPoints.size(),
    contactArea/scalar(refContactPoints.size())
);

if (useDirectionalContactQuadrature)
{
    refContactNormals =
        contactProperties.get<vectorField>("contactNormals");
    refContactPointAreas =
        contactProperties.get<scalarField>("contactPointAreas");

    if
    (
        refContactNormals.size() != refContactPoints.size()
     || refContactPointAreas.size() != refContactPoints.size()
    )
    {
        FatalErrorInFunction
            << "contactPoints, contactNormals and contactPointAreas must "
            << "have identical sizes" << exit(FatalError);
    }
}
const scalar foundationStiffness =
    contactProperties.get<scalar>("foundationStiffness");
const scalar contactDampingRatio =
    contactProperties.get<scalar>("dampingRatio");
const scalar staticFriction =
    contactProperties.get<scalar>("staticFriction");
const scalar dynamicFriction =
    contactProperties.get<scalar>("dynamicFriction");
const scalar frictionVelocity =
    contactProperties.get<scalar>("frictionVelocity");
const scalar bearingCapacity =
    contactProperties.get<scalar>("bearingCapacity");
const scalar contactSearchRadius =
    contactProperties.get<scalar>("searchRadius");
const scalar contactFluidDensity =
    contactProperties.get<scalar>("fluidDensity");

// A body-only retaining frame keeps a freely moving small specimen inside the
// short laboratory test section.  It does not alter the fluid boundaries.
const bool useBodyRetainer =
    immersedBodyProperties.found("retainerProperties");
vector retainerMin(-GREAT, -GREAT, -GREAT);
vector retainerMax(GREAT, GREAT, GREAT);
vector retainerActiveAxes(Zero);
scalar retainerStiffness = 0;
scalar retainerDampingRatio = 0;
scalar retainerMaxForce = GREAT;
if (useBodyRetainer)
{
    const dictionary& retainerProperties =
        immersedBodyProperties.subDict("retainerProperties");
    retainerMin = retainerProperties.get<vector>("min");
    retainerMax = retainerProperties.get<vector>("max");
    retainerActiveAxes = retainerProperties.get<vector>("activeAxes");
    retainerStiffness = retainerProperties.get<scalar>("stiffness");
    retainerDampingRatio =
        retainerProperties.get<scalar>("dampingRatio");
    retainerMaxForce = retainerProperties.get<scalar>("maxForce");
}

const dictionary& stabilityProperties =
    immersedBodyProperties.subDict("stabilityProperties");
const scalar kFieldMin = stabilityProperties.get<scalar>("kMin");
const scalar kFieldMax = stabilityProperties.get<scalar>("kMax");
const scalar omegaFieldMin = stabilityProperties.get<scalar>("omegaMin");
const scalar omegaFieldMax = stabilityProperties.get<scalar>("omegaMax");
const scalar nutFieldMax = stabilityProperties.get<scalar>("nutMax");

const dictionary& sixDoFCoeffs =
    immersedBodyProperties.subDict("sixDoFCoeffs");

sixDoFRigidBodyMotion bodyMotion
(
    sixDoFCoeffs,
    sixDoFCoeffs,
    runTime
);

volScalarField bodyMask
(
    IOobject
    (
        "bodyMask",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, Zero)
);

volVectorField bodyVelocity
(
    IOobject
    (
        "bodyVelocity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimVelocity, Zero)
);

vector hydroForce(Zero);
vector hydroMoment(Zero);
vector contactForce(Zero);
vector contactMoment(Zero);

fileName bodyMotionDirectory
(
    runTime.path()/"postProcessing"/"immersedBodyMotion"
);
autoPtr<OFstream> bodyMotionLogPtr;
if (Pstream::master())
{
    mkDir(bodyMotionDirectory);
    bodyMotionLogPtr.reset
    (
        new OFstream(bodyMotionDirectory/"motion.dat")
    );
    bodyMotionLogPtr()
        << "# time centreOfMass orientation velocity angularVelocity "
        << "hydroForce hydroMoment contactForce contactMoment" << nl;
}

#include "updateImmersedBodyFields.H"

// Avoid an artificial start-up impulse from water initially occupying the
// solid mask. The prescribed inlet flow remains unchanged outside the body.
U.primitiveFieldRef() *= scalar(1) - bodyMask.primitiveField();
U.correctBoundaryConditions();
