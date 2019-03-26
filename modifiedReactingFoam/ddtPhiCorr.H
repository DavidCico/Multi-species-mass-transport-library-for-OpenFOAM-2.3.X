#ifndef DDTPHICORR_H
#define DDTPHICORR_H


template<class Type>
tmp<GeometricField<typename flux<Type>::type, fvsPatchField, surfaceMesh> >
ddtCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField
    <
        typename flux<Type>::type,
        fvsPatchField,
        surfaceMesh
    >& phi
)
{
    return fv::ddtScheme<Type>::New
    (
        U.mesh(),
        U.mesh().ddtScheme("ddt(" + rho.name() + ',' + U.name() + ')')
    )().fvcDdtPhiCorr(rho, U, phi);
}




#endif // DDTPHICORR_H
