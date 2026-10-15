/**
 * SVG Hierarchy Tree and Wiring Form Inspector for Step 5 Bundle Viewer.
 */

export class BundleWiringManager {
    constructor(state, callbacks) {
        this.state = state;
        this.callbacks = callbacks || {};
    }

    renderGroupsTree(treeContainer, grpTree, activeGroupId, wiring) {
        if (!treeContainer) return;
        treeContainer.innerHTML = '';

        if (!grpTree || grpTree.length === 0) {
            treeContainer.innerHTML = `<div style="font-size: 0.72rem; color: #64748b; padding: 0.5rem;">No labeled groups found in SVG. Elements with inkscape:label will appear here.</div>`;
            return;
        }

        const renderNode = (node, depth) => {
            const item = document.createElement('div');
            item.className = `group-tree-item ${node.id === activeGroupId ? 'selected' : ''}`;
            item.style.paddingLeft = `${depth * 14 + 6}px`;

            const cfg = wiring[node.id] || wiring[node.label];
            let roleBadge = '';
            if (cfg) {
                if (cfg.role === 'arrow') roleBadge = '<span class="chip-mini blue">Arrow</span>';
                else if (cfg.role === 'state_molecule') roleBadge = '<span class="chip-mini purple">State+Mol</span>';
                else if (cfg.role === 'state') roleBadge = '<span class="chip-mini amber">State</span>';
                else if (cfg.role === 'molecule') roleBadge = '<span class="chip-mini emerald">Mol</span>';
                else if (cfg.role === 'push_arrow') roleBadge = '<span class="chip-mini purple">Push</span>';
            }

            item.innerHTML = `
                <div style="display: flex; align-items: center; gap: 4px; overflow: hidden; flex: 1;">
                    <span style="font-size: 0.7rem; color: #64748b;">${node.tag === 'g' ? '📁' : '🔹'}</span>
                    <span style="overflow: hidden; text-overflow: ellipsis; white-space: nowrap;">${node.label}</span>
                </div>
                ${roleBadge}
            `;

            item.onclick = (e) => {
                e.stopPropagation();
                if (this.callbacks.onSelectGroup) {
                    this.callbacks.onSelectGroup(node.id);
                }
            };
            treeContainer.appendChild(item);

            if (node.children && node.children.length > 0) {
                node.children.forEach(child => renderNode(child, depth + 1));
            }
        };

        grpTree.forEach(node => renderNode(node, 0));
    }

    renderWiringPanel(panel, grpId, cfg = {}) {
        if (!panel) return;
        const currentRole = cfg.role || 'none';
        const startFrame = cfg.startFrame !== undefined ? cfg.startFrame : 0;
        const endFrame = cfg.endFrame !== undefined ? cfg.endFrame : 100;
        const targetFrame = cfg.targetFrame !== undefined ? cfg.targetFrame : 0;
        const direction = cfg.direction || 'forward';
        const molName = cfg.moleculeName || grpId;

        panel.innerHTML = `
            <div style="font-weight: 700; color: #38bdf8; font-size: 0.8rem; margin-bottom: 0.35rem; display: flex; justify-content: space-between;">
                <span>Element: <code>${grpId}</code></span>
                <span class="chip blue" style="font-size: 0.68rem;">Role: ${currentRole.toUpperCase()}</span>
            </div>

            <!-- Role Selector Buttons -->
            <div style="display: flex; gap: 0.3rem; margin-bottom: 0.6rem; flex-wrap: wrap;">
                <button class="btn btn-xxs btn-role ${currentRole === 'none' ? 'active btn-primary' : 'btn-secondary'}" data-role="none">None</button>
                <button class="btn btn-xxs btn-role ${currentRole === 'arrow' ? 'active btn-primary' : 'btn-secondary'}" data-role="arrow">Reaction Arrow</button>
                <button class="btn btn-xxs btn-role ${currentRole === 'state' ? 'active btn-primary' : 'btn-secondary'}" data-role="state">State</button>
                <button class="btn btn-xxs btn-role ${currentRole === 'state_molecule' ? 'active btn-primary' : 'btn-secondary'}" data-role="state_molecule">State + Mol</button>
                <button class="btn btn-xxs btn-role ${currentRole === 'molecule' ? 'active btn-primary' : 'btn-secondary'}" data-role="molecule">Molecule</button>
                <button class="btn btn-xxs btn-role ${currentRole === 'push_arrow' ? 'active btn-primary' : 'btn-secondary'}" data-role="push_arrow">Push Arrow</button>
            </div>

            <!-- Dynamic Role Parameter Inputs -->
            <div id="role-fields-container" style="background: rgba(0,0,0,0.3); padding: 0.5rem; border-radius: 4px; border: 1px solid rgba(255,255,255,0.06); display: flex; flex-direction: column; gap: 0.4rem;">
                ${this.renderRoleFieldsHtml(currentRole, { startFrame, endFrame, targetFrame, direction, molName })}
            </div>

            <div style="display: flex; justify-content: flex-end; gap: 0.4rem; margin-top: 0.6rem;">
                <button id="btn-save-group-wiring" class="btn btn-xs btn-primary">Save Wiring</button>
            </div>
        `;

        const fieldsBox = panel.querySelector('#role-fields-container');
        panel.querySelectorAll('.btn-role').forEach(btn => {
            btn.addEventListener('click', () => {
                panel.querySelectorAll('.btn-role').forEach(b => b.classList.remove('active', 'btn-primary'));
                btn.classList.add('active', 'btn-primary');
                const newRole = btn.dataset.role;
                fieldsBox.innerHTML = this.renderRoleFieldsHtml(newRole, { startFrame, endFrame, targetFrame, direction, molName });
                this.bindRoleFieldEvents(panel, grpId);
            });
        });

        this.bindRoleFieldEvents(panel, grpId);

        panel.querySelector('#btn-save-group-wiring').addEventListener('click', () => {
            const roleBtn = panel.querySelector('.btn-role.active');
            const activeRole = roleBtn ? roleBtn.dataset.role : 'none';
            const newConfig = this.gatherWiringForm(panel, activeRole);
            if (this.callbacks.onSaveWiring) {
                this.callbacks.onSaveWiring(grpId, newConfig);
            }
        });
    }

    renderRoleFieldsHtml(role, data) {
        if (role === 'arrow') {
            return `
                <div style="display: flex; gap: 0.8rem; align-items: center; font-size: 0.74rem;">
                    <div>
                        <label>Start Frame:</label>
                        <input type="number" id="inp-arrow-start" value="${data.startFrame}" class="input-xs" style="width: 60px;">
                        <button class="btn btn-xxs btn-secondary" id="btn-use-cur-start">Cur</button>
                    </div>
                    <div>
                        <label>End Frame:</label>
                        <input type="number" id="inp-arrow-end" value="${data.endFrame}" class="input-xs" style="width: 60px;">
                        <button class="btn btn-xxs btn-secondary" id="btn-use-cur-end">Cur</button>
                    </div>
                    <div>
                        <label>Direction:</label>
                        <select id="sel-arrow-dir" class="input-xs">
                            <option value="forward" ${data.direction === 'forward' ? 'selected' : ''}>Forward ▶</option>
                            <option value="reverse" ${data.direction === 'reverse' ? 'selected' : ''}>Reverse ◀</option>
                        </select>
                    </div>
                </div>
            `;
        }
        if (role === 'state') {
            return `
                <div style="display: flex; gap: 0.8rem; align-items: center; font-size: 0.74rem;">
                    <div>
                        <label>Target Frame:</label>
                        <input type="number" id="inp-state-target" value="${data.targetFrame}" class="input-xs" style="width: 60px;">
                        <button class="btn btn-xxs btn-secondary" id="btn-use-cur-target">Cur</button>
                    </div>
                </div>
            `;
        }
        if (role === 'state_molecule') {
            return `
                <div style="display: flex; gap: 0.8rem; align-items: center; font-size: 0.74rem; flex-wrap: wrap;">
                    <div>
                        <label>Target Frame:</label>
                        <input type="number" id="inp-state-target" value="${data.targetFrame}" class="input-xs" style="width: 55px;">
                        <button class="btn btn-xxs btn-secondary" id="btn-use-cur-target">Cur</button>
                    </div>
                    <div>
                        <label>Molecule/Bundle:</label>
                        <input type="text" id="inp-mol-name" value="${data.molName}" class="input-xs" style="width: 110px;">
                    </div>
                    <button class="btn btn-xxs btn-secondary" id="btn-inspect-mol">Orbitals</button>
                </div>
            `;
        }
        if (role === 'molecule') {
            return `
                <div style="display: flex; gap: 0.8rem; align-items: center; font-size: 0.74rem;">
                    <div>
                        <label>Molecule Name:</label>
                        <input type="text" id="inp-mol-name" value="${data.molName}" class="input-xs" style="width: 140px;">
                    </div>
                    <button class="btn btn-xs btn-secondary" id="btn-inspect-mol">Open Orbitals</button>
                </div>
            `;
        }
        if (role === 'push_arrow') {
            return `
                <div style="display: flex; gap: 0.8rem; align-items: center; font-size: 0.74rem;">
                    <div>
                        <label>Active Start:</label>
                        <input type="number" id="inp-push-start" value="${data.startFrame}" class="input-xs" style="width: 60px;">
                    </div>
                    <div>
                        <label>Active End:</label>
                        <input type="number" id="inp-push-end" value="${data.endFrame}" class="input-xs" style="width: 60px;">
                    </div>
                    <div>
                        <label>Direction:</label>
                        <select id="sel-push-dir" class="input-xs">
                            <option value="forward" ${data.direction === 'forward' ? 'selected' : ''}>Forward Only ▶</option>
                            <option value="reverse" ${data.direction === 'reverse' ? 'selected' : ''}>Reverse Only ◀</option>
                            <option value="both" ${data.direction === 'both' ? 'selected' : ''}>Both Directions ◀▶</option>
                        </select>
                    </div>
                </div>
            `;
        }
        return `<div style="font-size: 0.72rem; color: #64748b;">No interactive action assigned.</div>`;
    }

    bindRoleFieldEvents(panel, grpId) {
        panel.querySelector('#btn-use-cur-start')?.addEventListener('click', () => {
            const inp = panel.querySelector('#inp-arrow-start');
            if (inp) inp.value = this.state.currentFrame;
        });
        panel.querySelector('#btn-use-cur-end')?.addEventListener('click', () => {
            const inp = panel.querySelector('#inp-arrow-end');
            if (inp) inp.value = this.state.currentFrame;
        });
        panel.querySelector('#btn-use-cur-target')?.addEventListener('click', () => {
            const inp = panel.querySelector('#inp-state-target');
            if (inp) inp.value = this.state.currentFrame;
        });
        panel.querySelector('#btn-inspect-mol')?.addEventListener('click', () => {
            const molName = panel.querySelector('#inp-mol-name')?.value || grpId;
            if (this.callbacks.onOpenMolecule) {
                this.callbacks.onOpenMolecule(grpId, { moleculeName: molName });
            }
        });
    }

    gatherWiringForm(panel, role) {
        if (role === 'none') return { role: 'none' };
        if (role === 'arrow') {
            return {
                role: 'arrow',
                startFrame: parseInt(panel.querySelector('#inp-arrow-start')?.value || 0, 10),
                endFrame: parseInt(panel.querySelector('#inp-arrow-end')?.value || 100, 10),
                direction: panel.querySelector('#sel-arrow-dir')?.value || 'forward'
            };
        }
        if (role === 'state') {
            return {
                role: 'state',
                targetFrame: parseInt(panel.querySelector('#inp-state-target')?.value || 0, 10)
            };
        }
        if (role === 'state_molecule') {
            return {
                role: 'state_molecule',
                targetFrame: parseInt(panel.querySelector('#inp-state-target')?.value || 0, 10),
                moleculeName: panel.querySelector('#inp-mol-name')?.value || 'Molecule'
            };
        }
        if (role === 'molecule') {
            return {
                role: 'molecule',
                moleculeName: panel.querySelector('#inp-mol-name')?.value || 'Molecule'
            };
        }
        if (role === 'push_arrow') {
            return {
                role: 'push_arrow',
                startFrame: parseInt(panel.querySelector('#inp-push-start')?.value || 0, 10),
                endFrame: parseInt(panel.querySelector('#inp-push-end')?.value || 100, 10),
                direction: panel.querySelector('#sel-push-dir')?.value || 'both'
            };
        }
        return { role: 'none' };
    }
}
