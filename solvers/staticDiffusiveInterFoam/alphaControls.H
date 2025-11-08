// Make the alpha1 dictionary accessible via the alphaControls keyword
const dictionary& alphaControls = mesh.solverDict(alpha1.name());

// Number of times to solve alphaEqn (updates on each loop)
label nAlphaCorr(alphaControls.lookup<label>("nAlphaCorr"));

// Alpha equation subcycling, useful where the interface time step is limiting 
// Must be equal to 1 - an if statement in alphaEqn.H checks this
label nAlphaSubCycles(alphaControls.lookup<label>("nAlphaSubCycles"));

// Semi-implicit MULES keyword (see OpenFOAM Foundation release notes 2.3.0)
// Must be false - an if statemnent in alphaEqn.H checks this 
bool MULESCorr(alphaControls.lookupOrDefault<Switch>("MULESCorr", false));
