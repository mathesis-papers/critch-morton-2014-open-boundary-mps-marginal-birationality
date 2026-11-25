"""Exact finite examples of coherent contraction and reconstruction."""

from fractions import Fraction
from itertools import product
import json
from pathlib import Path
import random

import verify


COVECTORS = ((1, 0), (0, 1), (1, 1), (2, -1))


def combination(matrices, covector):
    return tuple(tuple(sum(Fraction(covector[index]) * matrices[index][row][column]
                           for index in range(2))
                       for column in range(2)) for row in range(2))


def product_matrix(matrices, covectors):
    result = verify.IDENTITY
    for covector in covectors:
        result = verify.matrix_mul(result, combination(matrices, covector))
    return result


def covectors_for(length, start, retained):
    before = tuple(COVECTORS[index % len(COVECTORS)] for index in range(start))
    after = tuple(COVECTORS[(index + start + 1) % len(COVECTORS)]
                  for index in range(length - start - retained))
    return before, after


def contract(tensor, start, retained, before, after):
    output = dict.fromkeys(product(range(2), repeat=retained), Fraction(0))
    for word, value in tensor.items():
        coefficient = Fraction(1)
        for position, covector in enumerate(before):
            coefficient *= covector[word[position]]
        for position, covector in enumerate(after, start + retained):
            coefficient *= covector[word[position]]
        output[word[start:start + retained]] += coefficient * value
    return output


def weighted_reconstruct(length, data, before, after):
    left = verify.row_mul(data['left'], verify.inverse(product_matrix(data['matrices'], before)))
    right = verify.matrix_vec(verify.inverse(product_matrix(data['matrices'], after)), data['right'])
    return verify.state(length, data['matrices'], left, right)


def in_chart(tensor):
    kernel = tuple((tensor[first, middle, 0], tensor[first, middle, 1],
                    -tensor[0, first, middle], -tensor[1, first, middle])
                   for first, middle in product(range(2), repeat=2))
    minor = tuple(row[:3] for row in kernel[:3])
    divisor = verify.determinant(minor)
    if not divisor:
        return False
    cramer = tuple(verify.determinant(tuple(tuple(-kernel[row][3] if column == replaced
                                                else minor[row][column]
                                                for column in range(3)) for row in range(3)))
                   for replaced in range(3))
    raw_left = (cramer[2], divisor)
    basis = tuple(tuple(sum(raw_left[letter] * tensor[letter, row, column]
                            for letter in range(2)) for column in range(2)) for row in range(2))
    slices = tuple(tuple(tuple(tensor[row, letter, column] for column in range(2))
                          for row in range(2)) for letter in range(2))
    matrices = tuple(verify.matrix_mul(slice_matrix, verify.adjugate(basis))
                     for slice_matrix in slices)
    return bool(verify.determinant(basis) and verify.determinant(verify.matrix_add(*matrices))
                and verify.determinant(slices[0])
                and tensor[0, 0, 0] * tensor[0, 1, 1] - tensor[0, 0, 1] * tensor[0, 1, 0])


def admissible(matrices, left, right):
    if not in_chart(verify.state(3, matrices, left, right)):
        return False
    if not all(verify.determinant(combination(matrices, covector)) for covector in COVECTORS):
        return False
    for length in range(4, 8):
        tensor = verify.state(length, matrices, left, right)
        for start in range(length - 2):
            before, after = covectors_for(length, start, 3)
            if not in_chart(contract(tensor, start, 3, before, after)):
                return False
            if not in_chart(verify.marginal(tensor, start)):
                return False
        for retained in range(4, length + 1):
            for start in range(length - retained + 1):
                before, after = covectors_for(length, start, retained)
                block = contract(tensor, start, retained, before, after)
                if not in_chart(contract(block, 0, 3, (), ((1, 1),) * (retained - 3))):
                    return False
    return True


def check_sample(matrices, left, right):
    counters = {'triple_source_cases': 0, 'triple_target_cases': 0,
                'longer_block_cases': 0, 'coordinate_equalities': 0,
                'quartic_chart_cases': 0, 'seven_parameter_cases': 0}
    target = verify.state(3, matrices, left, right)
    target_data = verify.chart(target)
    for covector in COVECTORS:
        assert verify.determinant(combination(matrices, covector))
    for length in range(4, 8):
        tensor = verify.state(length, matrices, left, right)
        for start in range(length - 2):
            before, after = covectors_for(length, start, 3)
            local = contract(tensor, start, 3, before, after)
            data = verify.chart(local)
            recovered = weighted_reconstruct(length, data, before, after)
            assert recovered == tensor
            counters['triple_source_cases'] += 1
            counters['coordinate_equalities'] += len(tensor)
            extended = weighted_reconstruct(length, target_data, before, after)
            assert contract(extended, start, 3, before, after) == target
            counters['triple_target_cases'] += 1
            counters['coordinate_equalities'] += 8
            ordinary = verify.marginal(tensor, start)
            ordinary_data = verify.chart(ordinary)
            assert verify.polynomial_formula(length, start, ordinary_data) == tensor
            counters['quartic_chart_cases'] += 1
            counters['coordinate_equalities'] += len(tensor)
            normalized_left = data['left']
            assert normalized_left[1] == 1
            first, second = data['matrices']
            parameter = normalized_left[0]
            assert first[1] == (1 - parameter * first[0][0], -parameter * first[0][1])
            assert second[1] == (-parameter * second[0][0], 1 - parameter * second[0][1])
            assert verify.state(3, data['matrices'], normalized_left, data['right']) == local
            counters['seven_parameter_cases'] += 1
            counters['coordinate_equalities'] += 8
        for retained in range(4, length + 1):
            for start in range(length - retained + 1):
                before, after = covectors_for(length, start, retained)
                block = contract(tensor, start, retained, before, after)
                triple = contract(block, 0, 3, (), ((1, 1),) * (retained - 3))
                data = verify.chart(triple)
                combined_after = ((1, 1),) * (retained - 3) + after
                recovered = weighted_reconstruct(length, data, before, combined_after)
                assert recovered == tensor
                assert contract(recovered, start, retained, before, after) == block
                counters['longer_block_cases'] += 1
                counters['coordinate_equalities'] += len(tensor) + len(block)
    return counters


def main():
    generator = random.Random(20261103)
    candidates = [(((Fraction(1), Fraction(2)), (Fraction(3), Fraction(4))),
                   ((Fraction(5), Fraction(6)), (Fraction(7), Fraction(8))),
                   (Fraction(1), Fraction(0)), (Fraction(1), Fraction(2)))]
    accepted = []
    totals = None
    for attempt in range(100):
        if attempt == 0:
            first, second, left, right = candidates[0]
        else:
            first, second = tuple(tuple(tuple(Fraction(generator.randint(-4, 5))
                                             for column in range(2)) for row in range(2))
                                  for index in range(2))
            left, right = tuple(tuple(Fraction(generator.randint(-3, 4)) for index in range(2))
                                for boundary in range(2))
        if not admissible((first, second), left, right):
            continue
        counters = check_sample((first, second), left, right)
        accepted.append({'matrices': [[[str(entry) for entry in row] for row in matrix]
                                       for matrix in (first, second)],
                         'left': list(map(str, left)), 'right': list(map(str, right))})
        if totals is None:
            totals = dict.fromkeys(counters, 0)
        for key, count in counters.items():
            totals[key] += count
        if len(accepted) == 3:
            break
    assert len(accepted) == 3
    output = {'arithmetic': 'exact rational numbers; Python standard library',
              'lengths': [4, 5, 6, 7], 'parameter_samples': accepted,
              'covectors': COVECTORS, 'counts': totals,
              'scope': 'Finite examples on the stated rational charts; Corollaries 1-4 in note.tex have written proofs.',
              'status': 'all equalities hold'}
    text = json.dumps(output, indent=2) + '\n'
    Path(__file__).with_name('extension-checks.json').write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
