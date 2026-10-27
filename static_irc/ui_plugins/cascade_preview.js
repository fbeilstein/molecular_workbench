/**
 * Cascade Sequence List and Junction Inspector UI Builders.
 */

export function renderSequenceList(dom, state, onModified) {
    dom.seqList.innerHTML = '';
    const currentTrajs = state.currentTrajs;
    if (currentTrajs.length === 0) {
        const noSeq = document.createElement('div');
        noSeq.style.cssText = 'font-size:12px;color:#64748b;';
        noSeq.textContent = 'No previous trajectories found. Complete Step 4 first.';
        dom.seqList.appendChild(noSeq);
        return;
    }

    if (state.orderModified || currentTrajs.length !== (state.originalTrajs || []).length) {
        const notice = document.createElement('div');
        notice.style.cssText = 'font-size:11px;color:#d97706;background:#fef3c7;border:1px solid #fde68a;border-radius:4px;padding:4px 8px;margin-bottom:6px;display:flex;justify-content:space-between;align-items:center;';
        const noticeText = document.createElement('span');
        noticeText.textContent = 'Sequence modified. Click ⚡ Re-Stitch Pathway to apply.';
        notice.appendChild(noticeText);

        const resetBtn = document.createElement('button');
        resetBtn.textContent = '↺ Reset';
        resetBtn.style.cssText = 'padding:1px 6px;font-size:10px;border:1px solid #d97706;background:#fff;color:#b45309;border-radius:3px;cursor:pointer;';
        resetBtn.onclick = (e) => {
            e.stopPropagation();
            state.currentTrajs = [...(state.originalTrajs || [])];
            state.orderModified = false;
            renderSequenceList(dom, state, onModified);
            if (onModified) onModified();
        };
        notice.appendChild(resetBtn);
        dom.seqList.appendChild(notice);
    }

    currentTrajs.forEach((t, i) => {
        const sRow = document.createElement('div');
        sRow.style.cssText = 'display:flex;align-items:center;gap:6px;font-size:12px;padding:6px 8px;background:#ffffff;border:1px solid #e2e8f0;border-radius:6px;';

        const badge = document.createElement('span');
        badge.style.cssText = 'padding:2px 6px;border-radius:4px;font-size:11px;font-weight:600;background:#dbeafe;color:#1e40af;';
        badge.textContent = `Step ${i + 1}`;

        const sName = document.createElement('span');
        sName.style.cssText = 'flex:1;font-family:monospace;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;';
        const label = (t.reaction_smiles || t.filename || '').replace('>>', ' ➔ ');
        sName.title = label;
        sName.textContent = label;

        const fCount = document.createElement('span');
        fCount.style.cssText = 'color:#64748b;font-size:11px;';
        fCount.textContent = `${t.n_frames}f`;

        const btnBox = document.createElement('div');
        btnBox.style.cssText = 'display:flex;gap:2px;';

        const flipBtn = document.createElement('button');
        flipBtn.textContent = t.is_reversed ? '⇄ Inv' : '⇄';
        flipBtn.title = t.is_reversed ? 'Reaction is inverted' : 'Invert reaction direction';
        flipBtn.style.cssText = `padding:2px 5px;font-size:10px;line-height:1;border:${t.is_reversed ? '1px solid #f59e0b' : '1px solid #cbd5e1'};background:${t.is_reversed ? '#fef3c7' : '#f8fafc'};color:${t.is_reversed ? '#b45309' : '#334155'};border-radius:3px;cursor:pointer;`;
        flipBtn.onclick = (e) => {
            e.stopPropagation();
            t.is_reversed = !t.is_reversed;
            state.orderModified = true;
            renderSequenceList(dom, state, onModified);
            if (onModified) onModified();
        };

        const upBtn = document.createElement('button');
        upBtn.textContent = '▲';
        upBtn.style.cssText = `padding:2px 5px;font-size:10px;line-height:1;border:1px solid #cbd5e1;background:#f8fafc;border-radius:3px;cursor:${i > 0 ? 'pointer' : 'default'};opacity:${i > 0 ? '1' : '0.3'};`;
        upBtn.disabled = i === 0;
        upBtn.onclick = (e) => {
            e.stopPropagation();
            if (i > 0) {
                const temp = currentTrajs[i];
                currentTrajs[i] = currentTrajs[i - 1];
                currentTrajs[i - 1] = temp;
                state.orderModified = true;
                renderSequenceList(dom, state, onModified);
                if (onModified) onModified();
            }
        };

        const downBtn = document.createElement('button');
        downBtn.textContent = '▼';
        downBtn.style.cssText = `padding:2px 5px;font-size:10px;line-height:1;border:1px solid #cbd5e1;background:#f8fafc;border-radius:3px;cursor:${i < currentTrajs.length - 1 ? 'pointer' : 'default'};opacity:${i < currentTrajs.length - 1 ? '1' : '0.3'};`;
        downBtn.disabled = i === currentTrajs.length - 1;
        downBtn.onclick = (e) => {
            e.stopPropagation();
            if (i < currentTrajs.length - 1) {
                const temp = currentTrajs[i];
                currentTrajs[i] = currentTrajs[i + 1];
                currentTrajs[i + 1] = temp;
                state.orderModified = true;
                renderSequenceList(dom, state, onModified);
                if (onModified) onModified();
            }
        };

        const delBtn = document.createElement('button');
        delBtn.textContent = '✕';
        delBtn.style.cssText = 'padding:2px 5px;font-size:10px;line-height:1;border:1px solid #fca5a5;background:#fef2f2;color:#dc2626;border-radius:3px;cursor:pointer;';
        delBtn.onclick = (e) => {
            e.stopPropagation();
            currentTrajs.splice(i, 1);
            state.orderModified = true;
            renderSequenceList(dom, state, onModified);
            if (onModified) onModified();
        };

        btnBox.append(flipBtn, upBtn, downBtn, delBtn);
        sRow.append(badge, sName, fCount, btnBox);
        dom.seqList.appendChild(sRow);
    });
}

export function populateJunctions(dom, reports, junctionChoices) {
    dom.juncList.innerHTML = '';
    dom.juncSelect.innerHTML = '';

    if (!reports || reports.length === 0) {
        const noJ = document.createElement('div');
        noJ.style.cssText = 'font-size:12px;color:#94a3b8;';
        noJ.textContent = 'No intermediate junctions detected.';
        dom.juncList.appendChild(noJ);

        const opt = document.createElement('option');
        opt.textContent = 'No junctions';
        dom.juncSelect.appendChild(opt);
        return;
    }

    reports.forEach((rep, jIdx) => {
        const jBox = document.createElement('div');
        jBox.style.cssText = 'padding:8px;background:#ffffff;border:1px solid #cbd5e1;border-radius:6px;';

        const jHdr = document.createElement('div');
        jHdr.style.cssText = 'display:flex;justify-content:space-between;font-size:12px;font-weight:600;margin-bottom:4px;';
        jHdr.textContent = `Junction ${jIdx + 1}➔${jIdx + 2}: ${rep.hero_formula || 'Intermediate'}`;

        const sel = document.createElement('select');
        sel.style.cssText = 'width:100%;padding:5px 8px;font-size:11.5px;border-radius:4px;border:1px solid #cbd5e1;';

        (rep.all_candidates || []).forEach(cand => {
            const opt = document.createElement('option');
            opt.value = cand.id;
            opt.textContent = `${cand.label} (RMSD: ${cand.rmsd != null ? cand.rmsd.toFixed(4) : '--'} Å)`;
            if (junctionChoices[jIdx] === cand.id || (!junctionChoices[jIdx] && cand.id.includes('mol'))) {
                opt.selected = true;
                junctionChoices[jIdx] = cand.id;
            }
            sel.appendChild(opt);
        });

        sel.onchange = () => { junctionChoices[jIdx] = sel.value; };
        jBox.append(jHdr, sel);
        dom.juncList.appendChild(jBox);

        const jOpt = document.createElement('option');
        jOpt.value = jIdx.toString();
        jOpt.textContent = `Junction ${jIdx + 1}➔${jIdx + 2} (${rep.hero_formula || 'Intermediate'})`;
        dom.juncSelect.appendChild(jOpt);
    });
}
