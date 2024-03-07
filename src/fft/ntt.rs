use crate::field::field_element::FieldElement;
use crate::utils::bit_reverse_copy::bit_reverse_copy;

// source:
// paper - https://www.nayuki.io/page/number-theoretic-transform-integer-dft
// code  - https://www.nayuki.io/res/number-theoretic-transform-integer-dft/numbertheoretictransform.py
pub fn ntt<'a>(
    root: FieldElement<'a>,
    inputs: Vec<FieldElement<'a>>,
) -> Vec<FieldElement<'a>> {
    let field = inputs[0].field;

    let mut inputs = bit_reverse_copy(field.zero(), inputs);
    let n = inputs.len();

    let powtable = {
        let mut powtable = Vec::with_capacity(n / 2);
        let mut temp = field.one();
        for _ in 0..n / 2 {
            powtable.push(temp);
            temp = temp * root;
        }
        powtable
    };

    let mut size = 1;
    while size < n {
        size = size << 1;

        let halfsize = size / 2;
        let tablestep = n / size;

        for i in (0..n).step_by(size) {
            let mut k = 0;
            for j in i..i+halfsize {
                let l = j + halfsize;

                let even = inputs[j];
                let odd  = inputs[l] * powtable[k];
                inputs[j] = even + odd;
                inputs[l] = even - odd;

                k += tablestep;
            }
        }
    };

    inputs
}

pub fn intt<'a>(
    root: FieldElement<'a>,
    input: Vec<FieldElement<'a>>,
) -> Vec<FieldElement<'a>> {
    let orig_n = input.len();
    if orig_n < 2 {
        return input;
    }
    let n = orig_n.next_power_of_two();

    let field = input[0].field;
    let ninv = FieldElement::new(&field, n as u128).inverse();

    ntt(root.inverse(), input)
        .into_iter()
        .map(|tv| ninv * tv)
        .collect()
}

#[cfg(test)]
mod tests {
    use crate::fft::ntt::{intt, ntt};
    use crate::field::field::{Field, FIELD_PRIME};
    use crate::field::field_element::FieldElement;
    use crate::field::polynomial::Polynomial;

    #[test]
    fn test_ntt () {
        let field = Field::new(FIELD_PRIME);
        let n = 1 << 4;
        let primitive_root = field.primitive_nth_root(n);

        // values can be random
        let input = [10350860596407318609598574026175964133, 60692809610834653822383343680910625982, 223446197944610152228521360138425742723, 123599176902523769876954930401435714041, 233214499950980668770362073427851594143, 197530481770421435151547222505733630031, 6028204552208455457232478170590637777, 129106051215868132791440857107220454376, 46875137253396986423834480299002499296, 40573479539486208028801437611599580111, 177627388180112816822358878396956962568, 63754231379381382860231899477157171256, 213977912421556511151382836938765186268, 247295448209556494808801789962732329479, 198078312580458497833840274756537503682, 140348661180454074099943144751461445367]
            .into_iter()
            .map(|x| FieldElement::new(&field, x))
            .collect::<Vec<_>>();
        let values = [219013573292644897785762424283206192714, 28020178707455534238013018981848447223, 125720672179066355667363683873634014638, 9544075888957995047526079628773702483, 236214009288214032104373542256167121711, 203576991437594049347129945434757211067, 161303837601531457486430204397030363075, 8066037348193233635957451882263404827, 106698173671205255857026330656055139947, 205516443913240407551582667880265743260, 132452175458644240344865798387130681692, 14403130148933356826258737037147692544, 103258398926149853393925877736501914903, 241567637358481607032146874821122458208, 184640833807669488035490783312852642403, 79102880510147994351196921219409543952]
            .into_iter()
            .map(|x| FieldElement::new(&field, x))
            .collect::<Vec<_>>();
        assert_eq!(input.len() as u128, n);
        assert_eq!(values.len() as u128, n);

        let result = ntt(primitive_root, input.clone());
        assert_eq!(result, values);

        let poly = Polynomial::new(input);
        let result = poly.evaluate_domain(
            &(0..values.len())
                .map(|i| primitive_root ^ i)
                .collect::<Vec<_>>()
        );
        assert_eq!(result, values)
    }

    #[test]
    fn test_intt () {
        let field = Field::new(FIELD_PRIME);
        let n = 1 << 4;
        let primitive_root = field.primitive_nth_root(n);

        // values can be random
        let values = vec![159, 179, 197, 143, 198, 82, 100, 153, 45, 158, 154, 238, 46, 121, 148, 200]
            .into_iter()
            .map(|x| FieldElement::new(&field, x))
            .collect::<Vec<_>>();
        let coeffs = vec![2321, 46679697743149797158402415879589215379, 85767599764045409871854383990500128680, 170048455543476672374689900824216177289, 56517926799859326797837626323965682333, 150718635918560071455504820257610329093, 149093701728889244918633279335367822666, 266977550113122771518657412035427200127, 270497897142230380135924736767050120990, 63434915687244166391766073758524869310, 261359683971832165794823307314869483630, 172866549451408128829178953127270691728, 213979970342371053338087110443084438582, 83513730222590766426030683639871516493, 44774808819693939686538502893362807298, 127752053889369146389468687545690486361]
            .into_iter()
            .map(|x| FieldElement::new(&field, x))
            .collect::<Vec<_>>();
        assert_eq!(values.len() as u128, n);
        assert_eq!(coeffs.len() as u128, n);

        let result = ntt(primitive_root, values.clone());
        assert_eq!(coeffs, result);

        let result = intt(primitive_root, coeffs.clone());
        assert_eq!(values, result);
    }
}
