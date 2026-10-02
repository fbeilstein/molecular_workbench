window.WB = window.WB || {};

window.WB.initResizer = function() {
    const resizer = document.getElementById('dragMe');
    const leftSide = document.querySelector('.ketcher-panel');
    const rightSide = document.querySelector('.viewer-panel');

    if(!resizer || !leftSide || !rightSide) return;

    let x = 0;
    let leftWidth = 0;

    const mouseDownHandler = function(e) {
        x = e.clientX;
        leftWidth = leftSide.getBoundingClientRect().width;
        document.addEventListener('mousemove', mouseMoveHandler);
        document.addEventListener('mouseup', mouseUpHandler);
        
        // disable pointer events on iframes so they don't swallow mouse events
        const iframes = document.querySelectorAll('iframe');
        iframes.forEach(ifr => ifr.style.pointerEvents = 'none');
        document.getElementById('viewer-3d').style.pointerEvents = 'none';
    };

    const mouseMoveHandler = function(e) {
        const dx = e.clientX - x;
        const containerWidth = resizer.parentNode.getBoundingClientRect().width;
        const newLeftWidth = ((leftWidth + dx) * 100) / containerWidth;
        leftSide.style.flex = `0 0 ${newLeftWidth}%`;
        rightSide.style.flex = '1 1 0%';
    };

    const mouseUpHandler = function() {
        document.removeEventListener('mousemove', mouseMoveHandler);
        document.removeEventListener('mouseup', mouseUpHandler);
        
        // re-enable pointer events
        const iframes = document.querySelectorAll('iframe');
        iframes.forEach(ifr => ifr.style.pointerEvents = '');
        document.getElementById('viewer-3d').style.pointerEvents = '';
    };

    resizer.addEventListener('mousedown', mouseDownHandler);

    // Log Resizer
    const logResizer = document.getElementById('logDragMe');
    const logDrawer = document.getElementById('log-drawer');
    if (logResizer && logDrawer) {
        let startY = 0;
        let startHeight = 0;

        const logMouseDownHandler = function(e) {
            if (!logDrawer.classList.contains('open')) return;
            startY = e.clientY;
            startHeight = logDrawer.getBoundingClientRect().height;
            logDrawer.classList.add('dragging');
            document.addEventListener('mousemove', logMouseMoveHandler);
            document.addEventListener('mouseup', logMouseUpHandler);
            
            const iframes = document.querySelectorAll('iframe');
            iframes.forEach(ifr => ifr.style.pointerEvents = 'none');
            const viewer = document.getElementById('viewer-3d');
            if(viewer) viewer.style.pointerEvents = 'none';
        };

        const logMouseMoveHandler = function(e) {
            const dy = startY - e.clientY; // moving up increases height
            let newHeight = startHeight + dy;
            if (newHeight < 100) newHeight = 100;
            logDrawer.style.setProperty('--log-height', newHeight + 'px');
        };

        const logMouseUpHandler = function() {
            logDrawer.classList.remove('dragging');
            document.removeEventListener('mousemove', logMouseMoveHandler);
            document.removeEventListener('mouseup', logMouseUpHandler);
            
            const iframes = document.querySelectorAll('iframe');
            iframes.forEach(ifr => ifr.style.pointerEvents = '');
            const viewer = document.getElementById('viewer-3d');
            if(viewer) viewer.style.pointerEvents = '';
        };

        logResizer.addEventListener('mousedown', logMouseDownHandler);
    }
};
