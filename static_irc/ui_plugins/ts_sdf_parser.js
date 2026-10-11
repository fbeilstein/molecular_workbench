/**
 * Minimal SDF V2000 parser and position updater for 3D TS Editor.
 */

export function parseMolBlock(molblock) {
    if (!molblock) return { atoms: [], bonds: [] };
    const lines = molblock.split('\n');
    let countsIdx = lines.findIndex(l => l.includes('V2000'));
    if (countsIdx < 0) {
        if (lines.length < 4) return { atoms: [], bonds: [] };
        countsIdx = 3;
    }

    const countsLine = lines[countsIdx];
    const numAtoms = parseInt(countsLine.substring(0, 3).trim(), 10);
    const numBonds = parseInt(countsLine.substring(3, 6).trim(), 10);
    if (isNaN(numAtoms)) return { atoms: [], bonds: [] };

    const atoms = [];
    for (let i = 0; i < numAtoms; i++) {
        const line = lines[countsIdx + 1 + i];
        if (!line) continue;
        atoms.push({
            x: parseFloat(line.substring(0, 10).trim()),
            y: parseFloat(line.substring(10, 20).trim()),
            z: parseFloat(line.substring(20, 30).trim()),
            symbol: line.substring(31, 34).trim()
        });
    }

    const bonds = [];
    for (let i = 0; i < numBonds; i++) {
        const line = lines[countsIdx + 1 + numAtoms + i];
        if (!line) continue;
        bonds.push({
            from: parseInt(line.substring(0, 3).trim(), 10) - 1,
            to: parseInt(line.substring(3, 6).trim(), 10) - 1,
            type: parseInt(line.substring(6, 9).trim(), 10)
        });
    }

    return { atoms, bonds };
}

export function updateMolBlockPositions(molblock, group) {
    if (!molblock) return molblock;
    const lines = molblock.split('\n');
    let countsIdx = lines.findIndex(l => l.includes('V2000'));
    if (countsIdx < 0) {
        if (lines.length < 4) return molblock;
        countsIdx = 3;
    }

    const numAtoms = parseInt(lines[countsIdx].substring(0, 3).trim(), 10);
    if (isNaN(numAtoms)) return molblock;

    group.updateMatrixWorld(true);
    const worldMatrix = group.matrixWorld;
    const atomMeshes = group.userData.atomMeshes || [];

    const pad = (str, len) => str.length >= len ? str.substring(0, len) : ' '.repeat(len - str.length) + str;

    for (let i = 0; i < numAtoms; i++) {
        if (!atomMeshes[i]) continue;
        const worldPos = atomMeshes[i].position.clone().applyMatrix4(worldMatrix);
        const xs = pad(worldPos.x.toFixed(4), 10);
        const ys = pad(worldPos.y.toFixed(4), 10);
        const zs = pad(worldPos.z.toFixed(4), 10);
        const rest = lines[countsIdx + 1 + i].substring(30);
        lines[countsIdx + 1 + i] = xs + ys + zs + rest;
    }

    return lines.join('\n');
}
