import { dot, gramSchmidt, norm } from '../src/lattice';
import { lllReduce } from '../src/lll';
import { attackLWE, generateLWE } from '../src/bkz';

function assertTrue(condition: boolean, message: string): void {
  if (!condition) {
    throw new Error(message);
  }
}

function approxEqual(a: number, b: number, eps = 1e-9): boolean {
  return Math.abs(a - b) <= eps;
}

async function run(): Promise<void> {
  const gs = gramSchmidt([
    [3, 1],
    [1, 2],
  ]);
  assertTrue(approxEqual(gs.mu[1][0], 0.5), `Expected mu21=0.5, got ${gs.mu[1][0]}`);

  assertTrue(dot([1, 2], [3, 4]) === 11, 'dot product sanity check failed');
  assertTrue(approxEqual(norm([3, 4]), 5), 'norm sanity check failed');

  const lll = lllReduce([
    [19, 2],
    [7, 1],
  ]);
  const shortestNormSq = Math.min(...lll.reducedBasis.map((v) => dot(v, v)));
  assertTrue(shortestNormSq <= 5, `Expected short vector norm^2 <= 5, got ${shortestNormSq}`);

  const trials = 20;
  let exact = 0;
  for (let i = 0; i < trials; i += 1) {
    const inst = await generateLWE(4, 71, 2);
    const atk = attackLWE(inst.A, inst.b, 71, 4);
    if (atk.recovered && atk.recovered.join(',') === inst.secret.join(',')) {
      exact += 1;
    }
  }
  const successRate = exact / trials;
  assertTrue(successRate >= 0.8, `Expected toy attack success rate >= 0.8, got ${successRate}`);

  const hard = await generateLWE(20, 3329, 1);
  const hardAttack = attackLWE(hard.A, hard.b, 3329, 20);
  assertTrue(!hardAttack.success && hardAttack.recovered === null, 'Expected n=20 hard case to fail');

  console.log('verify:algo passed');
  console.log(`toy_success_rate=${successRate.toFixed(3)}`);
}

await run();
