
export const MappingUI = {
    render: (schemaContainer, ui, state, step_id) => {
        window.currentMapping = ui.initial_mapping ? Object.assign({}, ui.initial_mapping) : {};
            window.approvedMappings = {}; // Set of approved mapping keys
            let selectedNode = null;
            let hideApproved = false;
            
            // Buttons bar
            const btnBar = document.createElement('div');
            btnBar.style.marginBottom = '10px';
            const toggleBtn = document.createElement('button');
            toggleBtn.innerText = "Hide Approved Wires";
            toggleBtn.onclick = () => {
                hideApproved = !hideApproved;
                toggleBtn.innerText = hideApproved ? "Show Approved Wires" : "Hide Approved Wires";
                drawWires();
            };
            
            const approveAllBtn = document.createElement('button');
            approveAllBtn.innerText = "Approve All Remaining";
            approveAllBtn.style.marginLeft = '10px';
            approveAllBtn.onclick = () => {
                for (let k in window.currentMapping) {
                    window.approvedMappings[k] = true;
                }
                drawWires();
            };
            
            btnBar.appendChild(toggleBtn);
            btnBar.appendChild(approveAllBtn);
            schemaContainer.appendChild(btnBar);
            
            const container = document.createElement('div');
            container.style.position = 'relative';
            container.style.width = '100%';
            container.style.minHeight = '600px';
            container.style.overflow = 'hidden';
            container.style.border = '1px solid #ccc';
            container.style.backgroundColor = '#f4f4f4';
            
            const zoomContainer = document.createElement('div');
            zoomContainer.style.position = 'absolute';
            zoomContainer.style.top = '0';
            zoomContainer.style.left = '0';
            zoomContainer.style.width = '100%';
            zoomContainer.style.height = '100%';
            zoomContainer.style.transformOrigin = '0 0';
            
            let scale = 1;
            let panX = 0;
            let panY = 0;
            let isPanning = false;
            
            container.addEventListener('wheel', (e) => {
                e.preventDefault();
                const zoomFactor = e.deltaY > 0 ? 0.9 : 1.1;
                const rect = container.getBoundingClientRect();
                const mx = e.clientX - rect.left;
                const my = e.clientY - rect.top;
                
                panX = mx - (mx - panX) * zoomFactor;
                panY = my - (my - panY) * zoomFactor;
                scale *= zoomFactor;
                
                zoomContainer.style.transform = `translate(${panX}px, ${panY}px) scale(${scale})`;
            });
            
            container.addEventListener('mousedown', (e) => {
                if (e.target.tagName === 'circle' || e.target.tagName === 'path') return;
                isPanning = true;
                container.style.cursor = 'move';
            });
            
            window.addEventListener('mousemove', (e) => {
                if (isPanning) {
                    panX += e.movementX;
                    panY += e.movementY;
                    zoomContainer.style.transform = `translate(${panX}px, ${panY}px) scale(${scale})`;
                }
            });
            
            window.addEventListener('mouseup', () => { 
                if (isPanning) {
                    isPanning = false; 
                    container.style.cursor = 'default';
                }
            });
            
            // Layout is managed by direct positioning of wrappers now.
            
            const overlay = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
            overlay.style.position = 'absolute';
            overlay.style.top = '0';
            overlay.style.left = '0';
            overlay.style.width = '100%';
            overlay.style.height = '100%';
            overlay.style.overflow = 'visible';
            overlay.style.pointerEvents = 'none';
            overlay.style.zIndex = '100'; // Ensure wires render ABOVE molecule cards
            zoomContainer.appendChild(overlay);
            
            const buildSide = (side, data, startX) => {
                data.svgs.forEach((svgStr, fragIdx) => {
                    const wrapper = document.createElement('div');
                    wrapper.style.position = 'absolute';
                    wrapper.style.top = '0';
                    wrapper.style.left = '0';
                    wrapper.style.cursor = 'grab';
                    if (side === 'R') {
                        wrapper.style.backgroundColor = 'rgba(235, 245, 255, 0.9)'; // Light blue
                        wrapper.style.border = '2px solid #a3c2f0';
                    } else {
                        wrapper.style.backgroundColor = 'rgba(255, 240, 235, 0.9)'; // Light red/orange
                        wrapper.style.border = '2px solid #f0bda3';
                    }
                    wrapper.style.borderRadius = '8px';
                    wrapper.style.padding = '10px';
                    wrapper.style.boxShadow = '0 2px 8px rgba(0,0,0,0.2)';
                    wrapper.style.zIndex = '10'; // Keep cards below the overlay
                    wrapper.innerHTML = svgStr;
                    
                    let fx = startX;
                    let fy = 50 + fragIdx * 300; // stagger vertically to prevent overlap
                    wrapper.style.transform = `translate(${fx}px, ${fy}px)`;
                    let isDraggingFrag = false;
                    
                    wrapper.addEventListener('mousedown', (e) => {
                        if (e.target.tagName === 'circle') return;
                        e.stopPropagation();
                        isDraggingFrag = true;
                        wrapper.style.cursor = 'grabbing';
                    });
                    
                    window.addEventListener('mousemove', (e) => {
                        if (isDraggingFrag) {
                            fx += e.movementX / scale;
                            fy += e.movementY / scale;
                            wrapper.style.transform = `translate(${fx}px, ${fy}px)`;
                            drawWires();
                        }
                    });
                    
                    window.addEventListener('mouseup', () => {
                        if (isDraggingFrag) {
                            isDraggingFrag = false;
                            wrapper.style.cursor = 'grab';
                        }
                    });
                    
                    const svgEl = wrapper.querySelector('svg');
                    if (svgEl) {
                        const coords = data.coords[fragIdx];
                        for (let atomIdx in coords) {
                            const pt = coords[atomIdx];
                            const circle = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
                            circle.setAttribute('cx', pt.x);
                            circle.setAttribute('cy', pt.y);
                            circle.setAttribute('r', '15');
                            circle.setAttribute('fill', 'transparent');
                            circle.setAttribute('stroke', 'transparent');
                            circle.style.cursor = 'pointer';
                            circle.style.pointerEvents = 'all';
                            circle.id = `atom_${side}_${fragIdx}_${atomIdx}`;
                            circle.dataset.side = side;
                            circle.dataset.key = `${side}_${fragIdx}_${atomIdx}`;
                            
                            circle.addEventListener('mouseover', () => circle.setAttribute('stroke', 'blue'));
                            circle.addEventListener('mouseout', () => circle.setAttribute('stroke', 'transparent'));
                            
                            circle.addEventListener('click', (e) => {
                                const key = circle.dataset.key;
                                if (!selectedNode) {
                                    selectedNode = {key, el: circle};
                                    circle.setAttribute('fill', 'rgba(0,0,255,0.2)');
                                } else {
                                    // Make connection, remove conflicting
                                    const makeMap = (rk, pk) => {
                                        // Remove any existing mapping that uses rk or pk
                                        for (let k in window.currentMapping) {
                                            if (k === rk || window.currentMapping[k] === pk) {
                                                delete window.currentMapping[k];
                                                delete window.approvedMappings[k];
                                            }
                                        }
                                        window.currentMapping[rk] = pk;
                                    };
                                    
                                    if (selectedNode.key.startsWith('R') && key.startsWith('P')) {
                                        makeMap(selectedNode.key, key);
                                    } else if (selectedNode.key.startsWith('P') && key.startsWith('R')) {
                                        makeMap(key, selectedNode.key);
                                    }
                                    selectedNode.el.setAttribute('fill', 'transparent');
                                    selectedNode = null;
                                    drawWires();
                                }
                            });
                            svgEl.appendChild(circle);
                        }
                    }
                    zoomContainer.appendChild(wrapper);
                });
            };
            
            // Layout: Reactants on left (x=50), Products on right (x=600)
            buildSide('R', ui.reactants, 50);
            buildSide('P', ui.products, 600);
            
            container.appendChild(zoomContainer);
            schemaContainer.appendChild(container);
            
            const drawWires = () => {
                overlay.innerHTML = '';
                const zoomRect = zoomContainer.getBoundingClientRect();
                
                for (let rKey in window.currentMapping) {
                    const pKey = window.currentMapping[rKey];
                    const rEl = document.getElementById(`atom_${rKey}`);
                    const pEl = document.getElementById(`atom_${pKey}`);
                    
                    const isApproved = window.approvedMappings[rKey];
                    
                    if (isApproved) {
                        if (rEl) rEl.setAttribute('fill', 'rgba(0,255,0,0.1)');
                        if (pEl) pEl.setAttribute('fill', 'rgba(0,255,0,0.1)');
                        if (hideApproved) continue; // skip drawing wire entirely
                    } else {
                        if (rEl) rEl.setAttribute('fill', 'transparent');
                        if (pEl) pEl.setAttribute('fill', 'transparent');
                    }
                    
                    if (rEl && pEl) {
                        const rRect = rEl.getBoundingClientRect();
                        const pRect = pEl.getBoundingClientRect();
                        
                        // Inverse transform to local SVG coords
                        const x1 = (rRect.left + rRect.width/2 - zoomRect.left) / scale;
                        const y1 = (rRect.top + rRect.height/2 - zoomRect.top) / scale;
                        const x2 = (pRect.left + pRect.width/2 - zoomRect.left) / scale;
                        const y2 = (pRect.top + pRect.height/2 - zoomRect.top) / scale;
                        
                        const line = document.createElementNS('http://www.w3.org/2000/svg', 'path');
                        const cx = (x1 + x2) / 2;
                        line.setAttribute('d', `M ${x1} ${y1} C ${x1+100} ${y1} ${x2-100} ${y2} ${x2} ${y2}`);
                        
                        if (isApproved) {
                            line.setAttribute('stroke', 'rgba(0,180,0,0.5)'); // Green for approved
                        } else {
                            line.setAttribute('stroke', 'rgba(255,165,0,0.7)'); // Orange for unapproved
                        }
                        
                        line.setAttribute('stroke-width', '4');
                        line.setAttribute('fill', 'none');
                        line.style.pointerEvents = 'all';
                        line.style.cursor = 'pointer';
                        
                        line.addEventListener('mouseover', () => {
                            line.setAttribute('stroke-width', '6');
                            if (rEl) {
                                rEl.setAttribute('stroke', 'red');
                                rEl.setAttribute('stroke-width', '3');
                            }
                            if (pEl) {
                                pEl.setAttribute('stroke', 'red');
                                pEl.setAttribute('stroke-width', '3');
                            }
                        });
                        
                        line.addEventListener('mouseout', () => {
                            line.setAttribute('stroke-width', '4');
                            if (rEl) rEl.setAttribute('stroke', 'transparent');
                            if (pEl) pEl.setAttribute('stroke', 'transparent');
                        });
                        
                        line.addEventListener('click', (e) => {
                            if (e.shiftKey) {
                                delete window.currentMapping[rKey];
                                delete window.approvedMappings[rKey];
                            } else {
                                if (window.approvedMappings[rKey]) {
                                    delete window.approvedMappings[rKey];
                                } else {
                                    window.approvedMappings[rKey] = true;
                                }
                            }
                            drawWires();
                        });
                        
                        line.addEventListener('contextmenu', (e) => {
                            e.preventDefault();
                            delete window.currentMapping[rKey];
                            delete window.approvedMappings[rKey];
                            drawWires();
                        });
                        
                        overlay.appendChild(line);
                    }
                }
            };
            
            setTimeout(drawWires, 100);
            window.addEventListener('resize', drawWires);
        
    },
    
    getData: async () => {
        return { mapping: window.currentMapping };
    }
};
