// Load the GA script dynamically
var script = document.createElement('script');
script.src = 'https://www.googletagmanager.com/gtag/js?id=G-VQ50PWGDHN';
script.async = true;
document.head.appendChild(script);

// Initialize gtag
window.dataLayer = window.dataLayer || [];
function gtag(){dataLayer.push(arguments);}
gtag('js', new Date());
gtag('config', 'G-VQ50PWGDHN');
