Grid Refinement Team!!!

So far we have an interior residual calculation and mesh refinement implemented.

Next steps are:

    - Impliment edge residual calculations

    - Ensure that our code works as one module

        i.e. there should be one interface for the other groups to access our module
        
        it should give one output\\
        (we don't want there to be several awkward steps to refine the grid)

    - Ensure that our interface is compatible with the formats the other groups expect.
