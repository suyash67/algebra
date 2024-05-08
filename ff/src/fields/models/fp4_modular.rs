use super::quadratic_extension::{QuadExtConfig, QuadExtField};
use crate::{
    fields::{Fp2, Fp2Config},
    CyclotomicMultSubgroup, Zero,
};
use core::{marker::PhantomData, ops::Not};

pub trait Fp4ModularConfig: 'static + Send + Sync {
    type Fp2Config: Fp2Config;

    /// This *must* equal (0, 1);
    /// see [[DESD06, Section 5.1]](https://eprint.iacr.org/2006/471.pdf).
    const NONRESIDUE: Fp2<Self::Fp2Config>;

    /// Coefficients for the Frobenius automorphism.
    /// non_residue^((modulus^i-1)/4) for i=0,1,2,3
    const FROBENIUS_COEFF_FP4_C1: &'static [Fp2<Self::Fp2Config>];

    #[inline(always)]
    fn mul_fp2_by_nonresidue_in_place(fe: &mut Fp2<Self::Fp2Config>) -> &mut Fp2<Self::Fp2Config> {
        // see [[DESD06, Section 5.1]](https://eprint.iacr.org/2006/471.pdf).
        let new_c1 = fe.c0;
        Self::Fp2Config::mul_fp_by_nonresidue_in_place(&mut fe.c1);
        fe.c0 = fe.c1;
        fe.c1 = new_c1;
        fe
    }
}

pub struct Fp4ModularConfigWrapper<P: Fp4ModularConfig>(PhantomData<P>);

impl<P: Fp4ModularConfig> QuadExtConfig for Fp4ModularConfigWrapper<P> {
    type BasePrimeField = <P::Fp2Config as Fp2Config>::Fp;
    type BaseField = Fp2<P::Fp2Config>;
    type FrobCoeff = Fp2<P::Fp2Config>;

    const DEGREE_OVER_BASE_PRIME_FIELD: usize = 4;

    const NONRESIDUE: Self::BaseField = P::NONRESIDUE;

    const FROBENIUS_COEFF_C1: &'static [Self::BaseField] = P::FROBENIUS_COEFF_FP4_C1;

    #[inline(always)]
    fn mul_base_field_by_nonresidue_in_place(fe: &mut Self::BaseField) -> &mut Self::BaseField {
        P::mul_fp2_by_nonresidue_in_place(fe)
    }

    fn mul_base_field_by_frob_coeff(fe: &mut Self::BaseField, power: usize) {
        *fe *= Self::FROBENIUS_COEFF_C1[power % Self::DEGREE_OVER_BASE_PRIME_FIELD];
    }
}

pub type Fp4Modular<P> = QuadExtField<Fp4ModularConfigWrapper<P>>;

impl<P: Fp4ModularConfig> Fp4Modular<P> {
    pub fn mul_by_fp(&mut self, element: &<P::Fp2Config as Fp2Config>::Fp) {
        self.c0.mul_assign_by_fp(element);
        self.c1.mul_assign_by_fp(element);
    }

    pub fn mul_by_fp2(&mut self, element: &Fp2<P::Fp2Config>) {
        self.c0 *= element;
        self.c1 *= element;
    }
}

impl<P: Fp4ModularConfig> CyclotomicMultSubgroup for Fp4Modular<P> {
    const INVERSE_IS_FAST: bool = true;
    fn cyclotomic_inverse_in_place(&mut self) -> Option<&mut Self> {
        self.is_zero().not().then(|| {
            self.conjugate_in_place();
            self
        })
    }
}
