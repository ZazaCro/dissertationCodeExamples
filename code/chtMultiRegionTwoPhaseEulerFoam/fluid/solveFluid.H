if (finalIter)
{
    mesh.data::add("finalIteration", true);
}

if (frozenFlow)
{
    #include "EEqns.H"
}
else
{
    // Iteration index
    int iteration = 0;

    // Fluid loop
    while (iteration < nFluidLoop)
    {
        Info << "Fluid loop " << iteration + 1 << endl;

        fluid.solve();
        fluid.correct();

        #include "YEqns.H"

        if (faceMomentum)
        {
            #include "pUf/UEqns.H"
            #include "EEqns.H"
            #include "pUf/pEqn.H"
        }
        else
        {
            #include "pU/UEqns.H"
            #include "EEqns.H"
            #include "pU/pEqn.H"
        }

        fluid.correctKinematics();

        // Update alpha's for new U
        fluid.correctTurbulence();
        
        iteration++;
    }
}

if (finalIter)
{
    mesh.data::remove("finalIteration");
}
