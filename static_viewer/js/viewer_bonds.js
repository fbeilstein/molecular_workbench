/**
 * Chemical bond perception and MolBlock generation for standalone viewer.
 */

export const COV_RADII = {
  H: 0.31, C: 0.76, N: 0.71, O: 0.66, F: 0.57,
  P: 1.07, S: 1.05, Cl: 1.02, Br: 1.20, I: 1.39
};

export const MAX_VALENCE = { H: 1, F: 1, Cl: 1, Br: 1, I: 1, O: 3, N: 4, C: 4, S: 6, P: 5 };

export function computeBonds(symbols, coords) {
  const bonds = [];
  const n = Math.min(symbols.length, coords.length);
  const valences = new Array(n).fill(0);
  const candidates = [];
  const adj = Array.from({ length: n }, () => []);
  const distMat = {};

  for (let i = 0; i < n; i++) {
    const s1 = symbols[i];
    const r1 = COV_RADII[s1] || 0.75;
    for (let j = i + 1; j < n; j++) {
      const s2 = symbols[j];
      const r2 = COV_RADII[s2] || 0.75;
      if (s1 === 'H' && s2 === 'H' && n > 2) continue;

      const dx = coords[i][0] - coords[j][0];
      const dy = coords[i][1] - coords[j][1];
      const dz = coords[i][2] - coords[j][2];
      const d = Math.sqrt(dx * dx + dy * dy + dz * dz);
      const maxD = 1.28 * (r1 + r2);

      if (d >= 0.4 && d <= maxD) {
        adj[i].push(j);
        adj[j].push(i);
        distMat[`${i}_${j}`] = d;
        candidates.push({ i, j, s1, s2, d });
      }
    }
  }

  const getDist = (u, v) => distMat[u < v ? `${u}_${v}` : `${v}_${u}`] || 1.4;

  const ringDoubleBonds = new Set();
  const ringBondsAll = new Set();
  const cycles = [];
  const visited = new Set();
  const path = [];

  const findCycles = (v, start, depth) => {
    if (depth === 6) {
      if (adj[v].includes(start)) cycles.push([...path]);
      return;
    }
    for (const u of adj[v]) {
      if ((symbols[u] === 'C' || symbols[u] === 'N') && !visited.has(u)) {
        visited.add(u);
        path.push(u);
        findCycles(u, start, depth + 1);
        path.pop();
        visited.delete(u);
      }
    }
  };

  for (let i = 0; i < n; i++) {
    if (symbols[i] === 'C' || symbols[i] === 'N') {
      visited.add(i);
      path.push(i);
      findCycles(i, i, 1);
      path.pop();
      visited.delete(i);
    }
  }

  const uniqueRings = [];
  const seenRings = new Set();
  for (const c of cycles) {
    const sorted = [...c].sort((a, b) => a - b).join(',');
    if (!seenRings.has(sorted)) {
      seenRings.add(sorted);
      uniqueRings.push(c);
    }
  }

  for (const ring of uniqueRings) {
    let isConjugated = true;
    for (let k = 0; k < 6; k++) {
      if (getDist(ring[k], ring[(k + 1) % 6]) > 1.48) {
        isConjugated = false;
        break;
      }
    }
    if (!isConjugated) continue;

    const edges0 = [[ring[0], ring[1]], [ring[2], ring[3]], [ring[4], ring[5]]];
    const edges1 = [[ring[1], ring[2]], [ring[3], ring[4]], [ring[5], ring[0]]];
    const sum0 = edges0.reduce((acc, [u, v]) => acc + getDist(u, v), 0);
    const sum1 = edges1.reduce((acc, [u, v]) => acc + getDist(u, v), 0);
    const chosen = sum0 <= sum1 ? edges0 : edges1;

    for (const [u, v] of chosen) ringDoubleBonds.add(u < v ? `${u}_${v}` : `${v}_${u}`);
    for (let k = 0; k < 6; k++) ringBondsAll.add(ring[k] < ring[(k + 1) % 6] ? `${ring[k]}_${ring[(k + 1) % 6]}` : `${ring[(k + 1) % 6]}_${ring[k]}`);
  }

  for (const c of candidates) {
    const pairKey = c.i < c.j ? `${c.i}_${c.j}` : `${c.j}_${c.i}`;
    const pair = [c.s1, c.s2].sort().join('-');
    let targetOrder = 1;

    if (ringDoubleBonds.has(pairKey)) targetOrder = 2;
    else if (!ringBondsAll.has(pairKey)) {
      if (pair === 'C-O') targetOrder = c.d < 1.16 ? 3 : (c.d < 1.265 ? 2 : 1);
      else if (pair === 'C-C') targetOrder = c.d < 1.24 ? 3 : (c.d < 1.39 ? 2 : 1);
      else if (pair === 'C-N') targetOrder = c.d < 1.18 ? 3 : (c.d < 1.32 ? 2 : 1);
      else if (pair === 'N-N') targetOrder = c.d < 1.14 ? 3 : (c.d < 1.28 ? 2 : 1);
      else if (pair === 'N-O' || pair === 'O-O') targetOrder = c.d < 1.25 ? 2 : 1;
    }

    bonds.push({ a1: c.i, a2: c.j, order: 1, target: targetOrder, s1: c.s1, s2: c.s2, d: c.d, pairKey });
    valences[c.i]++;
    valences[c.j]++;
  }

  bonds.sort((a, b) => {
    const prioA = a.target === 3 ? 0 : (ringDoubleBonds.has(a.pairKey) ? 1 : (a.target === 2 ? 2 : 3));
    const prioB = b.target === 3 ? 0 : (ringDoubleBonds.has(b.pairKey) ? 1 : (b.target === 2 ? 2 : 3));
    return prioA !== prioB ? prioA - prioB : a.d - b.d;
  });

  for (const b of bonds) {
    if (b.target > 1) {
      const v1Cap = MAX_VALENCE[b.s1] || 4;
      const v2Cap = MAX_VALENCE[b.s2] || 4;
      const maxBoost = Math.min(v1Cap - valences[b.a1], v2Cap - valences[b.a2]);
      const boost = Math.max(0, Math.min(b.target - 1, maxBoost));
      if (boost > 0) {
        b.order += boost;
        valences[b.a1] += boost;
        valences[b.a2] += boost;
      }
    }
  }

  return bonds;
}

export function buildMolBlock(symbols, coords) {
  const natoms = Math.min(symbols.length, coords.length);
  const bonds = computeBonds(symbols, coords);
  let sdf = `\n  Reaction Pathway Viewer\n\n${String(natoms).padStart(3, ' ')}${String(bonds.length).padStart(3, ' ')}  0  0  0  0  0  0  0  0999 V2000\n`;
  for (let i = 0; i < natoms; i++) {
    const x = Number(coords[i][0]).toFixed(4).padStart(10, ' ');
    const y = Number(coords[i][1]).toFixed(4).padStart(10, ' ');
    const z = Number(coords[i][2]).toFixed(4).padStart(10, ' ');
    sdf += `${x}${y}${z} ${symbols[i].padEnd(3, ' ')} 0  0  0  0  0  0  0  0  0  0  0  0\n`;
  }
  for (let i = 0; i < bonds.length; i++) {
    sdf += `${String(bonds[i].a1 + 1).padStart(3, ' ')}${String(bonds[i].a2 + 1).padStart(3, ' ')}${String(bonds[i].order).padStart(3, ' ')}  0  0  0  0\n`;
  }
  sdf += 'M  END\n$$$$\n';
  return sdf;
}
