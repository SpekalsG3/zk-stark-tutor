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
                let left = inputs[j];
                let right = inputs[l] * powtable[k];
                inputs[j] = left + right;
                inputs[l] = left - right;
                k += tablestep;
            }
        }
    };

    inputs
}

#[cfg(test)]
mod tests {
    use crate::fft::ntt::ntt;
    use crate::field::field::{Field, FIELD_PRIME};
    use crate::field::field_element::FieldElement;
    use crate::field::polynomial::Polynomial;

    #[test]
    fn test () {
        let field = Field::new(FIELD_PRIME);
        let n = 1 << 4;
        let primitive_root = field.primitive_nth_root(n);

        let input = [10350860596407318609598574026175964133, 60692809610834653822383343680910625982, 223446197944610152228521360138425742723, 123599176902523769876954930401435714041, 233214499950980668770362073427851594143, 197530481770421435151547222505733630031, 6028204552208455457232478170590637777, 129106051215868132791440857107220454376, 46875137253396986423834480299002499296, 40573479539486208028801437611599580111, 177627388180112816822358878396956962568, 63754231379381382860231899477157171256, 213977912421556511151382836938765186268, 247295448209556494808801789962732329479, 198078312580458497833840274756537503682, 140348661180454074099943144751461445367]
            .into_iter()
            .map(|x| FieldElement::new(&field, x))
            .collect::<Vec<_>>();
        let output = [219013573292644897785762424283206192714, 28020178707455534238013018981848447223, 125720672179066355667363683873634014638, 9544075888957995047526079628773702483, 236214009288214032104373542256167121711, 203576991437594049347129945434757211067, 161303837601531457486430204397030363075, 8066037348193233635957451882263404827, 106698173671205255857026330656055139947, 205516443913240407551582667880265743260, 132452175458644240344865798387130681692, 14403130148933356826258737037147692544, 103258398926149853393925877736501914903, 241567637358481607032146874821122458208, 184640833807669488035490783312852642403, 79102880510147994351196921219409543952]
            .into_iter()
            .map(|x| FieldElement::new(&field, x))
            .collect::<Vec<_>>();

        let result = ntt(primitive_root, input.clone());
        assert_eq!(result, output);

        let poly = Polynomial::new(input);
        let result = poly.evaluate_domain(
            &(0..output.len())
                .map(|i| primitive_root ^ i)
                .collect::<Vec<_>>()
        );
        assert_eq!(output, result)
    }
}
