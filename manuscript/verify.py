"""Exact rational checks for the consecutive three-site reconstruction formula."""

from fractions import Fraction
from itertools import permutations, product
import json
from pathlib import Path
import random


IDENTITY = ((Fraction(1), Fraction(0)), (Fraction(0), Fraction(1)))
WORDS3 = tuple(product(range(2), repeat=3))


def matrix_add(first, second):
    return tuple(tuple(first[row][column] + second[row][column]
                       for column in range(2)) for row in range(2))


def matrix_mul(first, second):
    return tuple(tuple(sum(first[row][inner] * second[inner][column]
                           for inner in range(2))
                       for column in range(2)) for row in range(2))


def matrix_vec(matrix, vector):
    return tuple(sum(matrix[row][column] * vector[column] for column in range(2))
                 for row in range(2))


def row_mul(vector, matrix):
    return tuple(sum(vector[row] * matrix[row][column] for row in range(2))
                 for column in range(2))


def parity(order):
    return (-1) ** sum(order[first] > order[second]
                       for first in range(len(order))
                       for second in range(first + 1, len(order)))


def determinant(matrix):
    result = Fraction(0)
    for order in permutations(range(len(matrix))):
        term = Fraction(parity(order))
        for row, column in enumerate(order):
            term *= matrix[row][column]
        result += term
    return result


def adjugate(matrix):
    return ((matrix[1][1], -matrix[0][1]), (-matrix[1][0], matrix[0][0]))


def inverse(matrix):
    divisor = determinant(matrix)
    assert divisor
    return tuple(tuple(value / divisor for value in row) for row in adjugate(matrix))


def power(matrix, exponent):
    result = IDENTITY
    for _ in range(exponent):
        result = matrix_mul(result, matrix)
    return result


def amplitude(matrices, left, right, word):
    vector = left
    for letter in word:
        vector = row_mul(vector, matrices[letter])
    return sum(vector[index] * right[index] for index in range(2))


def state(length, matrices, left, right):
    return {word: amplitude(matrices, left, right, word)
            for word in product(range(2), repeat=length)}


def marginal(tensor, start):
    result = dict.fromkeys(WORDS3, Fraction(0))
    for word, value in tensor.items():
        result[word[start:start + 3]] += value
    return result


def chart(tensor):
    kernel = tuple((tensor[first, middle, 0], tensor[first, middle, 1],
                    -tensor[0, first, middle], -tensor[1, first, middle])
                   for first, middle in product(range(2), repeat=2))
    minor = tuple(row[:3] for row in kernel[:3])
    divisor = determinant(minor)
    assert divisor
    cramer = tuple(determinant(tuple(tuple(-kernel[row][3] if column == replaced
                                           else minor[row][column]
                                           for column in range(3)) for row in range(3)))
                   for replaced in range(3))
    raw_left = (cramer[2], divisor)
    raw_kernel = cramer[:2]
    raw_basis = tuple(tuple(sum(raw_left[letter] * tensor[letter, row, column]
                                for letter in range(2))
                            for column in range(2)) for row in range(2))
    slices = tuple(tuple(tuple(tensor[row, letter, column] for column in range(2))
                          for row in range(2)) for letter in range(2))
    raw_matrices = tuple(matrix_mul(slice_matrix, adjugate(raw_basis))
                         for slice_matrix in slices)
    raw_sum = matrix_add(*raw_matrices)
    raw_right = matrix_vec(raw_basis, raw_kernel)
    cross_minor = tensor[0, 0, 0] * tensor[0, 1, 1] - tensor[0, 0, 1] * tensor[0, 1, 0]
    factors = (divisor, determinant(raw_basis), determinant(raw_sum),
               determinant(slices[0]), cross_minor)
    assert all(factors)
    left = tuple(value / divisor for value in raw_left)
    basis = tuple(tuple(value / divisor for value in row) for row in raw_basis)
    matrices = tuple(matrix_mul(slice_matrix, inverse(basis)) for slice_matrix in slices)
    right = matrix_vec(basis, tuple(value / divisor for value in raw_kernel))
    for letter in range(2):
        assert row_mul(left, matrices[letter]) == tuple(Fraction(index == letter)
                                                      for index in range(2))
        assert matrix_vec(matrices[letter], right) == tuple(basis[index][letter]
                                                           for index in range(2))
    assert determinant(kernel) == 0
    kernel_vector = (*raw_kernel, *raw_left)
    assert all(sum(row[column] * kernel_vector[column] for column in range(4)) == 0
               for row in kernel)
    assert state(3, matrices, left, right) == tensor
    return {'matrices': matrices, 'left': left, 'right': right, 'factors': factors,
            'raw_matrices': raw_matrices, 'raw_sum': raw_sum, 'raw_left': raw_left,
            'raw_right': raw_right, 'cramer': cramer}


def reconstruct(length, start, data):
    sum_inverse = inverse(matrix_add(*data['matrices']))
    left = row_mul(data['left'], power(sum_inverse, start))
    right = matrix_vec(power(sum_inverse, length - 3 - start), data['right'])
    return state(length, data['matrices'], left, right)


def polynomial_formula(length, start, data):
    raw_adjugate = adjugate(data['raw_sum'])
    left = row_mul(data['raw_left'], power(raw_adjugate, start))
    right = matrix_vec(power(raw_adjugate, length - 3 - start), data['raw_right'])
    denominator = data['factors'][1] ** 3 * data['factors'][2] ** (length - 3)
    return {word: value / denominator for word, value in
            state(length, data['raw_matrices'], left, right).items()}


def quartic_coefficients():
    indexed_kernel = (
        ((1, 0), (1, 1), (-1, 0), (-1, 4)),
        ((1, 2), (1, 3), (-1, 1), (-1, 5)),
        ((1, 4), (1, 5), (-1, 2), (-1, 6)),
        ((1, 6), (1, 7), (-1, 3), (-1, 7)),
    )
    expansion = {}
    for order in permutations(range(4)):
        coefficient = parity(order)
        exponents = [0] * 8
        for row, column in enumerate(order):
            sign, index = indexed_kernel[row][column]
            coefficient *= sign
            exponents[index] += 1
        monomial = tuple(exponents)
        expansion[monomial] = expansion.get(monomial, 0) + coefficient
    expansion = {key: value for key, value in expansion.items() if value}
    source_terms = (
        (1, (3, 3, 4, 4)), (-1, (1, 3, 4, 5)), (-1, (2, 3, 4, 5)),
        (1, (0, 3, 5, 5)), (1, (1, 2, 3, 6)), (-1, (0, 3, 3, 6)),
        (-1, (2, 3, 4, 6)), (1, (1, 2, 5, 6)), (1, (1, 4, 5, 6)),
        (-1, (0, 5, 5, 6)), (-1, (1, 1, 6, 6)), (1, (0, 3, 6, 6)),
        (-1, (1, 2, 2, 7)), (1, (0, 2, 3, 7)), (1, (1, 1, 4, 7)),
        (1, (2, 2, 4, 7)), (-1, (0, 3, 4, 7)), (-1, (1, 4, 4, 7)),
        (-1, (0, 1, 5, 7)), (1, (0, 4, 5, 7)), (1, (0, 1, 6, 7)),
        (-1, (0, 2, 6, 7)),
    )
    source = {}
    for coefficient, indices in source_terms:
        exponents = tuple(indices.count(index) for index in range(8))
        source[exponents] = coefficient
    assert len(source) == len(expansion) == 22
    if expansion == source:
        return 1
    assert expansion == {key: -value for key, value in source.items()}
    return -1


def main():
    quartic_sign = quartic_coefficients()
    generator = random.Random(20140910)
    samples = [(((Fraction(1), Fraction(2)), (Fraction(3), Fraction(4))),
                ((Fraction(5), Fraction(6)), (Fraction(7), Fraction(8))),
                (Fraction(1), Fraction(0)), (Fraction(1), Fraction(2)))]
    attempts = 0
    while len(samples) < 6 and attempts < 100:
        attempts += 1
        matrices = tuple(tuple(tuple(Fraction(generator.randint(-4, 5))
                                     for _ in range(2)) for _ in range(2))
                         for _ in range(2))
        left = tuple(Fraction(generator.randint(-3, 4)) for _ in range(2))
        right = tuple(Fraction(generator.randint(-3, 4)) for _ in range(2))
        try:
            chart(state(3, matrices, left, right))
            for length in range(4, 9):
                tensor = state(length, matrices, left, right)
                for start in range(length - 2):
                    chart(marginal(tensor, start))
        except AssertionError:
            continue
        samples.append((*matrices, left, right))
    assert len(samples) == 6
    coordinate_equalities = 0
    source_cases = 0
    target_cases = 0
    for first, second, left, right in samples:
        matrices = (first, second)
        triple_data = chart(state(3, matrices, left, right))
        for length in range(4, 9):
            tensor = state(length, matrices, left, right)
            for start in range(length - 2):
                local = marginal(tensor, start)
                local_data = chart(local)
                recovered = reconstruct(length, start, local_data)
                assert recovered == tensor
                assert polynomial_formula(length, start, local_data) == recovered
                source_cases += 1
                coordinate_equalities += 2 * len(tensor)
                extended = reconstruct(length, start, triple_data)
                assert marginal(extended, start) == state(3, matrices, left, right)
                assert polynomial_formula(length, start, triple_data) == extended
                target_cases += 1
                coordinate_equalities += len(extended) + 8
    witness_tensor = state(3, samples[0][:2], samples[0][2], samples[0][3])
    assert tuple(witness_tensor[word] for word in WORDS3) == (145, 349, 337, 829, 489, 1173, 1193, 2933)
    witness = chart(witness_tensor)
    assert witness['factors'][0] == 387072
    assert witness['cramer'] == (118272, 10752, -1161216)
    output = {
        'arithmetic': 'exact integers and rational numbers; Python standard library',
        'quartic_coefficient_identity': {'terms': 22, 'determinant_equals_source_quartic_times': quartic_sign},
        'parameter_samples': 6,
        'lengths': [4, 5, 6, 7, 8],
        'windows': 'every consecutive triple at each listed length',
        'source_inverse_cases': source_cases,
        'target_inverse_cases': target_cases,
        'coordinate_equalities': coordinate_equalities,
        'witness_tensor_lexicographic': ['145', '349', '337', '829', '489', '1173', '1193', '2933'],
        'witness_chart_factors': [str(value) for value in witness['factors']],
        'witness_cramer_numerators': [str(value) for value in witness['cramer']],
        'scope': 'The quartic coefficient identity is exhaustive; reconstruction evaluations cover only the listed samples and lengths. The general theorem has a written proof in note.tex.',
        'status': 'all equalities hold',
    }
    text = json.dumps(output, indent=2) + '\n'
    Path(__file__).with_name('checks.json').write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
