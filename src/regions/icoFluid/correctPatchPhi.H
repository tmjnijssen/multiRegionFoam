        if (closedVolume_)
        {
            label intPatchID_ =
                mesh().boundaryMesh().findPatchID("interfaceShadow");

            forAll(phi_().boundaryField(), patchI)
            {
                if
                (
                    !phi_().boundaryField()[patchI].coupled()
                 && patchI != intPatchID_
//                 && isA<zeroGradientFvPatchVectorField>
//                    (
//                        U_().boundaryField()[patchI]
//                    )
                )
                {
                    phi_().boundaryField()[patchI] ==
                    (
                        U_().boundaryField()[patchI]
                        & mesh().Sf().boundaryField()[patchI]
                    );
                }
            }
        }

        if (!closedVolume_)
        {
            forAll(phi_().boundaryField(), patchI)
            {
                if
                (
                   !phi_().boundaryField()[patchI].coupled()
                 && !p_().boundaryField()[patchI].fixesValue()
                )
                {
                    phi_().boundaryField()[patchI] ==
                    (
                        U_().boundaryField()[patchI]
                        & mesh().Sf().boundaryField()[patchI]
                    );
                }
            }
        }
