"""CSS stylesheet for standalone bundle viewer export."""

STANDALONE_CSS = """
        :root {
            --bg-body: #0b0f17;
            --bg-surface: #111827;
            --border-color: rgba(255, 255, 255, 0.08);
            --primary: #38bdf8;
            --text-main: #f8fafc;
            --text-muted: #94a3b8;
        }
        * { box-sizing: border-box; margin: 0; padding: 0; }
        body {
            font-family: 'Inter', -apple-system, sans-serif;
            background: var(--bg-body);
            color: var(--text-main);
            min-height: 100vh;
            display: flex;
            flex-direction: column;
            padding: 1rem 1.5rem;
            gap: 1rem;
        }
        .header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            background: rgba(17, 24, 39, 0.85);
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 0.75rem 1.25rem;
        }
        .chips { display: flex; gap: 0.5rem; margin-top: 0.25rem; }
        .chip {
            background: rgba(255,255,255,0.06);
            border: 1px solid var(--border-color);
            padding: 2px 8px;
            border-radius: 9999px;
            font-size: 0.75rem;
        }
        .chip.blue { color: #38bdf8; }
        .chip.amber { color: #f59e0b; }
        .chip.emerald { color: #10b981; }
        .transport-bar {
            display: flex;
            align-items: center;
            gap: 0.6rem;
            background: rgba(0,0,0,0.3);
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 0.5rem 0.8rem;
        }
        .btn {
            background: rgba(255,255,255,0.08);
            border: 1px solid var(--border-color);
            color: #fff;
            padding: 6px 12px;
            border-radius: 6px;
            cursor: pointer;
            font-size: 0.82rem;
            font-weight: 600;
            transition: all 0.15s ease;
        }
        .btn:hover { background: rgba(255,255,255,0.16); }
        .btn-primary { background: #2563eb; border-color: #3b82f6; }
        .btn-primary:hover { background: #1d4ed8; }
        .btn-amber { background: #d97706; border-color: #f59e0b; }
        .dir-badge {
            font-size: 0.75rem;
            padding: 3px 8px;
            border-radius: 4px;
            font-weight: 700;
        }
        .dir-fwd { background: rgba(56,189,248,0.15); color: #38bdf8; border: 1px solid rgba(56,189,248,0.3); }
        .dir-rev { background: rgba(245,158,11,0.15); color: #f59e0b; border: 1px solid rgba(245,158,11,0.3); }
        .workspace {
            display: flex;
            gap: 1rem;
            height: 480px;
        }
        .viewer-box {
            flex: 1;
            background: var(--bg-surface);
            border: 1px solid var(--border-color);
            border-radius: 8px;
            position: relative;
            overflow: hidden;
        }
        .mechanism-box {
            flex: 1.1;
            background: var(--bg-surface);
            border: 1px solid var(--border-color);
            border-radius: 8px;
            display: flex;
            flex-direction: column;
            overflow: hidden;
        }
        .mechanism-header {
            padding: 0.6rem 1rem;
            border-bottom: 1px solid var(--border-color);
            font-size: 0.85rem;
            font-weight: 600;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }
        .svg-container {
            flex: 1;
            display: flex;
            align-items: center;
            justify-content: center;
            padding: 0.5rem;
            background: #ffffff;
            overflow: auto;
        }
        .svg-container svg {
            max-width: 100%;
            max-height: 100%;
        }
        .chart-box {
            height: 140px;
            background: var(--bg-surface);
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 0.75rem;
            position: relative;
        }
        .clickable-arrow {
            cursor: pointer;
            filter: drop-shadow(0 0 3px rgba(56,189,248,0.4));
            transition: all 0.2s ease;
        }
        .clickable-arrow:hover {
            filter: drop-shadow(0 0 8px rgba(56,189,248,0.9));
        }
        .clickable-mol {
            cursor: pointer;
        }
        .clickable-mol:hover rect, .clickable-mol:hover path {
            filter: drop-shadow(0 0 5px rgba(245,158,11,0.8));
        }
"""
