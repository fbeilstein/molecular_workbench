"""JavaScript playback and interactive logic for standalone bundle viewer export."""

STANDALONE_JS = """
        const MEGDATA = __JSON_DATA__;
        const WIRING = __WIRING_DATA__;
        const TOTAL_FRAMES = MEGDATA.frames ? MEGDATA.frames.length : 0;
        let curFrame = 0;
        let playDirection = 'forward';
        let isPlaying = false;
        let timer = null;

        // Init 3Dmol viewer
        const viewerEl = document.getElementById('viewer3d');
        const viewer = $3Dmol.createViewer(viewerEl, { backgroundColor: '#090d16' });
        let hasCamInit = false;

        function renderFrame(f) {
            curFrame = f;
            document.getElementById('scrubber').value = f;
            document.getElementById('readout').textContent = `Frame ${f} / ${TOTAL_FRAMES - 1}`;
            
            // Render 3D Molecule
            const frame = MEGDATA.frames[f];
            if (frame) {
                viewer.clear();
                let n = MEGDATA.symbols.length;
                let xyz = n + '\\nFrame ' + f + '\\n';
                for (let i = 0; i < n; i++) {
                    const s = MEGDATA.symbols[i];
                    const [x, y, z] = frame.coords[i];
                    xyz += s + ' ' + x.toFixed(4) + ' ' + y.toFixed(4) + ' ' + z.toFixed(4) + '\\n';
                }
                viewer.addModel(xyz, 'xyz');
                viewer.setStyle({}, {
                    stick: { radius: 0.14, colorscheme: 'Jmol' },
                    sphere: { scale: 0.28, colorscheme: 'Jmol' }
                });

                // Display selected active highres orbital if present
                if (frame.orbitals && frame.orbitals.active_cube) {
                    viewer.addVolumetricData(frame.orbitals.active_cube, "cube", {
                        isoval: 0.035,
                        color: "#38bdf8",
                        opacity: 0.8
                    });
                    viewer.addVolumetricData(frame.orbitals.active_cube, "cube", {
                        isoval: -0.035,
                        color: "#ef4444",
                        opacity: 0.8
                    });
                }

                if (!hasCamInit) {
                    viewer.zoomTo();
                    hasCamInit = true;
                }
                viewer.render();
            }

            drawPes();
        }

        // Wire interactive SVG arrows and molecules
        function wireInteractiveSvg() {
            const svgEl = document.querySelector('#svg-view svg');
            if (!svgEl) return;

            // Wire mapped reaction arrows
            if (WIRING.reactions) {
                Object.entries(WIRING.reactions).forEach(([rxnId, r]) => {
                    const el = svgEl.getElementById(rxnId) || svgEl.querySelector(`[data-id="${rxnId}"]`);
                    if (el) {
                        el.classList.add('clickable-arrow');
                        el.title = `Click to play: ${r.start_frame} → ${r.end_frame} (Forward)\nRight-click for Reverse`;
                        el.onclick = (e) => {
                            e.preventDefault();
                            playSegment(r.start_frame, r.end_frame, 'forward');
                        };
                        el.oncontextmenu = (e) => {
                            e.preventDefault();
                            playSegment(r.start_frame, r.end_frame, 'reverse');
                        };
                    }
                });
            }

            // Wire mapped molecules
            if (WIRING.molecules) {
                Object.entries(WIRING.molecules).forEach(([molId, m]) => {
                    const el = svgEl.getElementById(molId) || svgEl.querySelector(`[data-id="${molId}"]`);
                    if (el) {
                        el.classList.add('clickable-mol');
                        el.title = `${m.name || molId} (Frame ${m.frame}) - Click to jump`;
                        el.onclick = (e) => {
                            e.preventDefault();
                            pause();
                            renderFrame(m.frame);
                        };
                    }
                });
            }
        }

        function setDirection(dir) {
            playDirection = dir;
            const badge = document.getElementById('dir-badge');
            if (dir === 'forward') {
                badge.className = 'dir-badge dir-fwd';
                badge.textContent = 'Forward ⏩';
            } else {
                badge.className = 'dir-badge dir-rev';
                badge.textContent = 'Reverse ⏪';
            }
        }

        function playForward() {
            pause();
            setDirection('forward');
            isPlaying = true;
            timer = setInterval(() => {
                if (curFrame >= TOTAL_FRAMES - 1) {
                    pause();
                } else {
                    renderFrame(curFrame + 1);
                }
            }, 80);
        }

        function playReverse() {
            pause();
            setDirection('reverse');
            isPlaying = true;
            timer = setInterval(() => {
                if (curFrame <= 0) {
                    pause();
                } else {
                    renderFrame(curFrame - 1);
                }
            }, 80);
        }

        function pause() {
            isPlaying = false;
            if (timer) { clearInterval(timer); timer = null; }
        }

        function playSegment(sf, ef, dir) {
            pause();
            setDirection(dir);
            curFrame = (dir === 'forward') ? sf : ef;
            renderFrame(curFrame);
            isPlaying = true;
            timer = setInterval(() => {
                if (dir === 'forward') {
                    if (curFrame >= ef) pause();
                    else renderFrame(curFrame + 1);
                } else {
                    if (curFrame <= sf) pause();
                    else renderFrame(curFrame - 1);
                }
            }, 90);
        }

        // Draw PES Canvas
        const pesCanvas = document.getElementById('pes-canvas');
        const pesCtx = pesCanvas.getContext('2d');
        function drawPes() {
            const w = pesCanvas.parentElement.clientWidth;
            const h = pesCanvas.parentElement.clientHeight;
            if (pesCanvas.width !== w) pesCanvas.width = w;
            if (pesCanvas.height !== h) pesCanvas.height = h;

            pesCtx.clearRect(0, 0, w, h);
            const ep = (MEGDATA.energy_profile && MEGDATA.energy_profile.frame_energies_kcal) || [];
            if (!ep.length) return;

            const minE = Math.min(...ep) - 1.0;
            const maxE = Math.max(...ep) + 2.0;
            const padX = 35, padY = 20;
            const pw = w - padX * 2, ph = h - padY * 2;

            // Draw line
            pesCtx.beginPath();
            pesCtx.strokeStyle = '#38bdf8';
            pesCtx.lineWidth = 2.5;
            ep.forEach((e, i) => {
                const x = padX + (i / (ep.length - 1)) * pw;
                const y = h - padY - ((e - minE) / (maxE - minE)) * ph;
                if (i === 0) pesCtx.moveTo(x, y); else pesCtx.lineTo(x, y);
            });
            pesCtx.stroke();

            // Cursor
            const cx = padX + (curFrame / (ep.length - 1)) * pw;
            const cy = h - padY - ((ep[curFrame] - minE) / (maxE - minE)) * ph;
            pesCtx.beginPath();
            pesCtx.arc(cx, cy, 6, 0, Math.PI * 2);
            pesCtx.fillStyle = '#ef4444';
            pesCtx.fill();
            pesCtx.strokeStyle = '#fff';
            pesCtx.lineWidth = 2;
            pesCtx.stroke();
        }

        document.getElementById('btn-play-fwd').onclick = playForward;
        document.getElementById('btn-play-rev').onclick = playReverse;
        document.getElementById('btn-pause').onclick = pause;
        document.getElementById('scrubber').oninput = (e) => {
            pause();
            renderFrame(+e.target.value);
        };

        wireInteractiveSvg();
        renderFrame(0);
        window.addEventListener('resize', drawPes);
"""
