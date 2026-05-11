(function () {
    const prefetched = new Set();

    function prefetchUrl(url) {
        if (!url || prefetched.has(url)) return;
        prefetched.add(url);

        const link = document.createElement('link');
        link.rel = 'prefetch';
        link.href = url;
        document.head.appendChild(link);
    }

    function waitForDestination(url, timeoutMs) {
        return new Promise((resolve) => {
            let settled = false;
            const finish = () => {
                if (settled) return;
                settled = true;
                resolve();
            };

            const timer = window.setTimeout(finish, timeoutMs);
            const iframe = document.createElement('iframe');
            iframe.tabIndex = -1;
            iframe.setAttribute('aria-hidden', 'true');
            iframe.style.cssText = 'position:fixed;width:1px;height:1px;opacity:0;pointer-events:none;left:-10px;top:-10px;';
            iframe.onload = () => {
                window.clearTimeout(timer);
                iframe.remove();
                finish();
            };
            iframe.onerror = () => {
                window.clearTimeout(timer);
                iframe.remove();
                finish();
            };
            iframe.src = url;
            document.body.appendChild(iframe);
        });
    }

    async function prefetchAndNavigate(url) {
        prefetchUrl(url);
        const overlay = document.getElementById('white-transition');
        if (overlay) overlay.classList.add('is-active');

        await new Promise((resolve) => window.setTimeout(resolve, 430));
        await waitForDestination(url, 3000);
        window.location.href = url;
    }

    window.ZigBioPreloader = {
        prefetchUrl,
        prefetchAndNavigate,
    };
})();
