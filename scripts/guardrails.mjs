import { readFileSync, readdirSync, statSync } from 'node:fs';
import { join } from 'node:path';

function walk(dir) {
  const out = [];
  for (const entry of readdirSync(dir)) {
    const path = join(dir, entry);
    const st = statSync(path);
    if (st.isDirectory()) {
      if (entry === 'node_modules' || entry === 'dist' || entry === '.git') {
        continue;
      }
      out.push(...walk(path));
    } else {
      out.push(path);
    }
  }
  return out;
}

function fail(message) {
  console.error(message);
  process.exit(1);
}

const root = process.cwd();
const files = walk(root).filter((p) => {
  if (!/\.(ts|js|mjs|css|html|md)$/i.test(p)) {
    return false;
  }
  return !p.endsWith('scripts/guardrails.mjs');
});

for (const file of files) {
  const text = readFileSync(file, 'utf8');
  if (text.includes('Math.random(')) {
    fail(`Guardrail failed: Math.random found in ${file}`);
  }
}

const readme = readFileSync(join(root, 'README.md'), 'utf8').toLowerCase();
if (readme.includes('lll breaks kyber') || readme.includes('breaks nist')) {
  fail('Guardrail failed: prohibited claim found in README.md');
}

console.log('verify:guardrails passed');
