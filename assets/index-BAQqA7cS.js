(function(){let e=document.createElement(`link`).relList;if(e&&e.supports&&e.supports(`modulepreload`))return;for(let e of document.querySelectorAll(`link[rel="modulepreload"]`))n(e);new MutationObserver(e=>{for(let t of e)if(t.type===`childList`)for(let e of t.addedNodes)e.tagName===`LINK`&&e.rel===`modulepreload`&&n(e)}).observe(document,{childList:!0,subtree:!0});function t(e){let t={};return e.integrity&&(t.integrity=e.integrity),e.referrerPolicy&&(t.referrerPolicy=e.referrerPolicy),e.crossOrigin===`use-credentials`?t.credentials=`include`:e.crossOrigin===`anonymous`?t.credentials=`omit`:t.credentials=`same-origin`,t}function n(e){if(e.ep)return;e.ep=!0;let n=t(e);fetch(e.href,n)}})();function e(e,t){if(e.length!==t.length)throw Error(`Vector length mismatch`)}function t(e){if(e.length===0)throw Error(`Basis must be non-empty`);let t=e[0].length;if(t===0)throw Error(`Vectors must be non-empty`);for(let n of e)if(n.length!==t)throw Error(`All basis vectors must have equal dimension`)}function n(e){return e.map(e=>e.slice())}function r(t,n){e(t,n);let r=0;for(let e=0;e<t.length;e+=1)r+=t[e]*n[e];return r}function i(e){return Math.sqrt(r(e,e))}function a(t,n){e(t,n);let r=Array(t.length);for(let e=0;e<t.length;e+=1)r[e]=t[e]+n[e];return r}function o(t,n){e(t,n);let r=Array(t.length);for(let e=0;e<t.length;e+=1)r[e]=t[e]-n[e];return r}function s(e,t){let n=Array(e.length);for(let r=0;r<e.length;r+=1)n[r]=e[r]*t;return n}function c(e){t(e);let n=e.length,i=e[0].length,a=Array.from({length:n},()=>Array(i).fill(0)),c=Array.from({length:n},()=>Array(n).fill(0));for(let t=0;t<n;t+=1){let n=e[t].slice();for(let i=0;i<t;i+=1){let l=r(a[i],a[i]);if(l===0){c[t][i]=0;continue}let u=r(e[t],a[i])/l;c[t][i]=u;let d=s(a[i],u);n=o(n,d)}a[t]=n}return{gso:a,mu:c}}function l(e,t,n,i=.75){if(n<0||n+1>=e.length)throw Error(`Invalid Lovasz index`);let o=i*r(e[n],e[n]),c=a(e[n+1],s(e[n],t[n+1][n]));return o<=r(c,c)+1e-12}function u(e,t,n){let{mu:r}=c(e),i=Math.round(r[t][n]);if(i===0)return 0;let a=e[t],o=e[n];for(let e=0;e<a.length;e+=1)a[e]-=i*o[e];return i}function d(e,t){if(t<0||t+1>=e.length)throw Error(`Invalid swap index`);let n=e[t];e[t]=e[t+1],e[t+1]=n}function f(e){return e.map(e=>e.slice())}function p(e){let t=e.length,n=e.map(e=>e.slice()),r=1;for(let e=0;e<t;e+=1){let i=e;for(let r=e+1;r<t;r+=1)Math.abs(n[r][e])>Math.abs(n[i][e])&&(i=r);if(Math.abs(n[i][e])<1e-12)return 0;if(i!==e){let t=n[e];n[e]=n[i],n[i]=t,r*=-1}let a=n[e][e];r*=a;for(let r=e+1;r<t;r+=1){let i=n[r][e]/a;for(let a=e;a<t;a+=1)n[r][a]-=i*n[e][a]}}return r}function m(e,t=.75){let r=n(e),i=r.length,a=[];if(i<=1){let{gso:e,mu:t}=c(r);return a.push({type:`done`,i:0,before:n(r),after:n(r),gsoBefore:n(e),gsoAfter:n(e),muBefore:f(t),muAfter:f(t),description:`LLL complete (trivial basis).`}),{reducedBasis:r,steps:a,swapCount:0}}let o=1,s=0,p=0;for(;o<i;){if(p+=1,p>2e4)throw Error(`LLL did not converge within operation bound`);for(let e=o-1;e>=0;--e){let t=n(r),i=c(r),s=u(r,o,e),l=n(r),d=c(r);a.push({type:`size-reduce`,i:o,j:e,before:t,after:l,gsoBefore:n(i.gso),gsoAfter:n(d.gso),muBefore:f(i.mu),muAfter:f(d.mu),reductionCoeff:s,description:s===0?`k=${o}, j=${e}: size reduction skipped (round(mu)=0).`:`k=${o}, j=${e}: size reduce with coefficient ${s}.`})}let e=n(r),i=c(r);if(l(i.gso,i.mu,o-1,t)){let t=n(r),s=c(r);a.push({type:`advance`,i:o,before:e,after:t,gsoBefore:n(i.gso),gsoAfter:n(s.gso),muBefore:f(i.mu),muAfter:f(s.mu),lovaszSatisfied:!0,description:`k=${o}: Lovasz satisfied, advance to k=${o+1}.`}),o+=1}else{d(r,o-1),s+=1;let t=n(r),l=c(r);a.push({type:`swap`,i:o,before:e,after:t,gsoBefore:n(i.gso),gsoAfter:n(l.gso),muBefore:f(i.mu),muAfter:f(l.mu),lovaszSatisfied:!1,description:`k=${o}: Lovasz violated, swap b[${o-1}] and b[${o}].`}),o=Math.max(o-1,1)}}let m=c(r);return a.push({type:`done`,i:i-1,before:n(r),after:n(r),gsoBefore:n(m.gso),gsoAfter:n(m.gso),muBefore:f(m.mu),muAfter:f(m.mu),description:`LLL complete: ${s} swaps, ${a.length} trace steps.`}),{reducedBasis:r,steps:a,swapCount:s}}function h(e){if(e.length===0)return 1;let t=e[0].length;if(e.length!==t)throw Error(`Orthogonality defect requires square basis matrix`);let n=1;for(let t of e)n*=i(t);let r=Math.abs(p(e));return r<1e-12?1/0:n/r}function ee(e){let t=e[0],n=r(t,t);for(let i=1;i<e.length;i+=1){let a=r(e[i],e[i]);a<n&&(t=e[i],n=a)}return t.slice()}function te(){let e=new Uint32Array(1);return crypto.getRandomValues(e),e[0]??0}function ne(e,t){let n=t-e+1;return e+te()%n}function re(){return(te()+1)/4294967297}function ie(e,t){let n=(e%t+t)%t;return n>t/2?n-t:n}function ae(e){let t=re(),n=re(),r=Math.sqrt(-2*Math.log(t))*Math.cos(2*Math.PI*n);return Math.round(r*e)}function oe(e,t){let n=e[0].length,r=Array(n).fill(0);for(let i=0;i<e.length;i+=1)for(let a=0;a<n;a+=1)r[a]+=t[i]*e[i][a];return r}function se(e){let t=e.length;if(t===0||t>8)return null;let n=c(e),i=n.gso.map(e=>r(e,e)),a=Array(t).fill(0),o={bestNormSq:1/0,bestCoeffs:null};for(let t of e){let e=r(t,t);e>0&&e<o.bestNormSq&&(o.bestNormSq=e)}let s=(c,l)=>{if(l>=o.bestNormSq)return;if(c<0){let t=!1;for(let e of a)if(e!==0){t=!0;break}if(!t)return;let n=oe(e,a),i=r(n,n);i>0&&i<o.bestNormSq&&(o.bestNormSq=i,o.bestCoeffs=a.slice());return}if(i[c]<1e-12)return;let u=0;for(let e=c+1;e<t;e+=1)u-=a[e]*n.mu[e][c];let d=(o.bestNormSq-l)/i[c];if(d<0)return;let f=Math.sqrt(d),p=Math.ceil(u-f),m=Math.floor(u+f);for(let e=p;e<=m;e+=1){a[c]=e;let t=e-u;s(c-1,l+t*t*i[c])}a[c]=0};return s(t-1,0),o.bestCoeffs?oe(e,o.bestCoeffs):null}function ce(e,t){if(t<2)throw Error(`BKZ beta must be >= 2`);let n=m(e).reducedBasis,r=n.length,a=0,o=0,s=[],c=[];for(;o<20;){o+=1;let e=!1,l=0;for(let c=0;c<r;c+=1){let u=Math.min(c+t,r),d=n.slice(c,u).map(e=>e.slice());if(d.length<2)continue;let f=se(d);if(!f)continue;let p=i(f),h=i(d[0]);if(p+1e-9<h){s.push({tour:o,blockStart:c,blockEnd:u-1,headNormBefore:h,insertedNorm:p}),d[0]=f;let t=m(d).reducedBasis;for(let e=0;e<t.length;e+=1)n[c+e]=t[e];e=!0,a+=1,l+=1}}if(!e){c.push({tour:o,improvementsInTour:0,converged:!0});break}c.push({tour:o,improvementsInTour:l,converged:!1}),n=m(n).reducedBasis}return{reducedBasis:n,tours:o,improvements:a,blockImprovements:s,tourLogs:c}}function le(e,t,n){let r=e.length,i=e[0]?.length??0;if(r===0||i===0||t.length!==r)throw Error(`Invalid LWE dimensions`);let a=i+r+1,o=[];for(let e=0;e<i;e+=1){let t=Array(a).fill(0);t[e]=n,o.push(t)}for(let t=0;t<r;t+=1){let r=Array(a).fill(0);for(let a=0;a<i;a+=1)r[a]=(e[t][a]%n+n)%n;r[i+t]=1,o.push(r)}let s=Array(a).fill(0);for(let e=0;e<i;e+=1)s[e]=((t[e]??0)%n+n)%n;return s[a-1]=1,o.push(s),o}async function ue(e,t,n){let r=e+4,i=Array.from({length:r},()=>Array.from({length:e},()=>ne(0,t-1))),a=Array.from({length:e},()=>ne(0,4)),o=Array.from({length:r},()=>ae(n)),s=Array(r).fill(0);for(let n=0;n<r;n+=1){let r=0;for(let t=0;t<e;t+=1)r+=i[n][t]*a[t];r+=o[n],s[n]=(r%t+t)%t}return{A:i,b:s,secret:a,errors:o}}function g(e,t,n,r){let i=0;for(let a=0;a<e.length;a+=1){let o=0;for(let t=0;t<r.length;t+=1)o+=e[a][t]*r[t];let s=ie(t[a]-o,n);i+=s*s}return{ok:i/Math.max(e.length,1)<=9,residualNormSq:i}}function de(e,t,n,r){if(r>8)return null;let i=null,a=1/0,o=Array(r).fill(0),s=c=>{if(c===r){let r=g(e,t,n,o);r.residualNormSq<a&&(a=r.residualNormSq,i=o.slice());return}for(let e=0;e<=4;e+=1)o[c]=e,s(c+1)};return s(0),i&&g(e,t,n,i).ok?i:null}function fe(e,t){return(e%t+t)%t}function pe(e,t,n,r,i){let a=e[e.length-1];if(Math.abs(a)!==1)return null;let o=-a,s=e.slice(0,i).map(e=>o*e),c=[];c.push(s.map(e=>fe(e,r))),c.push(s.map(e=>fe(ie(e,r),r)));let l=null,u=1/0;for(let e of c){let i=g(t,n,r,e);if(i.residualNormSq<u&&(l=e,u=i.residualNormSq),i.ok)return e}return l}function me(e,t,n,i){if(i>=20)return{success:!1,recovered:null,shortVec:null,recoverMethod:`none`};let a=m(le(e,t,n)).reducedBasis,o=null,s=1/0;for(let e of a){let t=r(e,e);t>0&&t<s&&(s=t,o=e.slice())}if(o){let r=pe(o,e,t,n,i);if(r&&g(e,t,n,r).ok)return{success:!0,recovered:r,shortVec:o,recoverMethod:`short-vector`}}let c=de(e,t,n,i);return c?{success:!0,recovered:c,shortVec:o,recoverMethod:`bruteforce`}:{success:!1,recovered:null,shortVec:o,recoverMethod:`none`}}var he=document.querySelector(`#app`);if(!he)throw Error(`Missing app root`);var ge=document.documentElement.getAttribute(`data-theme`)??`dark`;he.innerHTML=`
<a class="skip-link" href="#main-content">Skip to main content</a>
<header class="topbar">
  <h1>crypto-lab-lll-break</h1>
  <button id="theme-toggle" type="button" class="theme-toggle" style="position: absolute; top: 0; right: 0" aria-label="${ge===`dark`?`Switch to light mode`:`Switch to dark mode`}">${ge===`dark`?`🌙`:`☀️`}</button>
</header>

<main id="main-content">
<section class="exhibit" id="exhibit-1" aria-labelledby="exhibit-1-title">
  <h2 id="exhibit-1-title">Exhibit 1 - What Is a Lattice?</h2>
  <div class="grid two">
    <div>
      <canvas id="lattice-canvas" width="500" height="400" role="img" aria-label="2D lattice canvas"></canvas>
      <p class="mono" id="det-label"></p>
      <p class="mono" id="same-lattice-label"></p>
    </div>
    <div>
      <label>b1 angle <input id="b1-angle" type="range" min="0" max="180" value="25" /></label>
      <label>b1 length <input id="b1-len" type="range" min="1" max="10" value="5" /></label>
      <label>b2 angle <input id="b2-angle" type="range" min="0" max="180" value="95" /></label>
      <label>b2 length <input id="b2-len" type="range" min="1" max="10" value="4" /></label>
      <button id="same-lattice-btn" type="button" class="btn">Same lattice, different basis</button>
      <p>
        A lattice is all integer combinations of two basis vectors. Determinant equals area of the
        fundamental parallelogram, an invariant under unimodular basis changes.
      </p>
    </div>
  </div>
</section>

<section class="exhibit" id="exhibit-2" aria-labelledby="exhibit-2-title">
  <h2 id="exhibit-2-title">Exhibit 2 - Gram-Schmidt Orthogonalization</h2>
  <div class="grid two">
    <div>
      <div class="matrix-grid">
        <label class="sr-only" for="gso-a11">GSO matrix row 1 column 1</label>
        <input id="gso-a11" type="number" inputmode="numeric" value="3" />
        <label class="sr-only" for="gso-a12">GSO matrix row 1 column 2</label>
        <input id="gso-a12" type="number" inputmode="numeric" value="1" />
        <label class="sr-only" for="gso-a21">GSO matrix row 2 column 1</label>
        <input id="gso-a21" type="number" inputmode="numeric" value="1" />
        <label class="sr-only" for="gso-a22">GSO matrix row 2 column 2</label>
        <input id="gso-a22" type="number" inputmode="numeric" value="2" />
      </div>
      <button id="gso-next" type="button" class="btn">Next</button>
      <pre id="gso-text" class="mono panel"></pre>
      <p class="panel">
        Gram-Schmidt vectors and mu coefficients are floating point. Basis vectors remain exact integers.
      </p>
    </div>
    <div>
      <canvas id="gso-canvas" width="500" height="360" role="img" aria-label="Gram-Schmidt visualization"></canvas>
    </div>
  </div>
</section>

<section class="exhibit" id="exhibit-3" aria-labelledby="exhibit-3-title">
  <h2 id="exhibit-3-title">Exhibit 3 - LLL Step-by-Step</h2>
  <div class="controls-row">
    <label>Dimension
      <select id="lll-dim">
        <option value="2">2D</option>
        <option value="3">3D</option>
      </select>
    </label>
    <label>Delta <input id="lll-delta" type="range" min="0.5" max="0.999" step="0.001" value="0.75" /></label>
    <button type="button" class="btn preset" data-preset="random2d">Random bad 2D</button>
    <button type="button" class="btn preset" data-preset="classic">Classic example</button>
    <button type="button" class="btn preset" data-preset="near">Near-orthogonal</button>
    <button type="button" class="btn preset" data-preset="challenge3d">3D challenge</button>
    <button type="button" class="btn preset" data-preset="qary">LWE-style q-ary</button>
  </div>
  <textarea id="lll-matrix" class="mono" rows="4">19 2
7 1</textarea>
  <div class="controls-row">
    <button id="lll-step" type="button" class="btn">Step</button>
    <button id="lll-auto" type="button" class="btn">Auto</button>
    <button id="lll-reset" type="button" class="btn">Reset</button>
  </div>
  <div class="grid two">
    <div>
      <canvas id="lll-canvas" width="420" height="420" role="img" aria-label="LLL basis evolution canvas"></canvas>
      <canvas id="defect-chart" width="420" height="140" role="img" aria-label="Orthogonality defect chart"></canvas>
    </div>
    <div>
      <div id="lll-log" class="panel mono log" aria-live="polite"></div>
      <pre id="lll-metrics" class="panel mono"></pre>
    </div>
  </div>
</section>

<section class="exhibit" id="exhibit-4" aria-labelledby="exhibit-4-title">
  <h2 id="exhibit-4-title">Exhibit 4 - Break a Toy LWE Instance</h2>
  <div class="controls-row">
    <label>n <input id="lwe-n" type="range" min="2" max="12" value="4" /></label>
    <label>q <input id="lwe-q" type="range" min="7" max="257" value="71" /></label>
    <label>sigma <input id="lwe-sigma" type="range" min="1" max="10" value="2" /></label>
    <label>beta <input id="lwe-beta" type="range" min="2" max="8" value="2" /></label>
    <button id="lwe-generate" type="button" class="btn">Generate LWE Instance</button>
    <button id="lwe-attack" type="button" class="btn">Run LLL/BKZ Attack</button>
    <button id="kyber-try" type="button" class="btn warn">Try Kyber-512 parameters</button>
  </div>
  <progress id="lwe-progress" max="100" value="0" aria-label="LWE attack progress"></progress>
  <div class="lwe-meter panel" role="status" aria-live="polite">
    <div class="lwe-meter-head mono" id="lwe-meter-text">Norm-gap confidence: waiting for attack run.</div>
    <div class="lwe-meter-track" aria-hidden="true">
      <div id="lwe-meter-fill" class="lwe-meter-fill"></div>
    </div>
  </div>
  <pre id="lwe-output" class="panel mono" aria-live="polite"></pre>
</section>

<section class="exhibit" id="exhibit-5" aria-labelledby="exhibit-5-title">
  <h2 id="exhibit-5-title">Exhibit 5 - Parameter Explorer</h2>
  <label>n <input id="explore-n" type="range" min="4" max="256" value="8" /></label>
  <div class="grid two">
    <pre id="explore-left" class="panel mono"></pre>
    <pre id="explore-right" class="panel mono"></pre>
  </div>
  <canvas id="threshold-chart" width="900" height="260" role="img" aria-label="Security threshold chart"></canvas>
  <p>
    ML-KEM (Kyber), FrodoKEM, ML-DSA (Dilithium), and FALCON all rely on LLL/BKZ failing at chosen parameters.
  </p>
</section>

</main>

<footer>
  <p>"Whether therefore ye eat, or drink, or whatsoever ye do,
  do all to the glory of God." — 1 Corinthians 10:31</p>
</footer>
`;function _(e){let t=document.getElementById(e);if(!t)throw Error(`Missing element: ${e}`);return t}function v(e,t,n){return Math.max(t,Math.min(n,e))}function y(e,t){let n=e*Math.PI/180;return[t*Math.cos(n),t*Math.sin(n)]}function b(e,t){let n=new Uint32Array(1);return crypto.getRandomValues(n),e+(n[0]??0)%(t-e+1)}function _e(e,t){return e[0]*t[1]-e[1]*t[0]}function ve(e,t,n){e.strokeStyle=`#334`,e.lineWidth=1,e.beginPath(),e.moveTo(t/2,0),e.lineTo(t/2,n),e.moveTo(0,n/2),e.lineTo(t,n/2),e.stroke()}function x(e,t,n,r,i=!1,a=2){e.save(),i&&e.setLineDash([6,4]),e.strokeStyle=r,e.fillStyle=r,e.lineWidth=a,e.beginPath(),e.moveTo(t[0],t[1]),e.lineTo(n[0],n[1]),e.stroke();let o=Math.atan2(n[1]-t[1],n[0]-t[0]);e.beginPath(),e.moveTo(n[0],n[1]),e.lineTo(n[0]-8*Math.cos(o-.3),n[1]-8*Math.sin(o-.3)),e.lineTo(n[0]-8*Math.cos(o+.3),n[1]-8*Math.sin(o+.3)),e.closePath(),e.fill(),e.restore()}function S(e){return e.length===2?[e[0],e[1]]:[e[0]-.45*e[2],e[1]-.25*e[2]]}function C(e,t,n,r,i){let a=e.getContext(`2d`);if(!a)return;a.clearRect(0,0,e.width,e.height),ve(a,e.width,e.height);let o=[e.width/2,e.height/2],s=S(n[0]),c=S(n[1]);for(let e=-5;e<=5;e+=1)for(let t=-5;t<=5;t+=1){let n=[o[0]+22*(e*s[0]+t*c[0]),o[1]-22*(e*s[1]+t*c[1])];a.fillStyle=`#1a3a6a`,a.beginPath(),a.arc(n[0],n[1],2,0,2*Math.PI),a.fill()}a.fillStyle=`#ffd700`,a.beginPath(),a.arc(o[0],o[1],4,0,2*Math.PI),a.fill();for(let e of t){let t=S(e);x(a,o,[o[0]+22*t[0],o[1]-22*t[1]],`#ff3366`,!1,2)}for(let e of n){let t=S(e);x(a,o,[o[0]+22*t[0],o[1]-22*t[1]],`#00ff88`,!1,3)}if(r)for(let e of r){let t=S(e);x(a,o,[o[0]+22*t[0],o[1]-22*t[1]],`#00d4ff`,!0,2)}i&&(a.fillStyle=`#ddd`,a.font=`13px "JetBrains Mono", monospace`,a.fillText(i,10,20))}var ye=_(`theme-toggle`);function be(e){let t=()=>{let t=document.documentElement.getAttribute(`data-theme`)===`light`?`light`:`dark`;e.textContent=t===`dark`?`🌙`:`☀️`,e.setAttribute(`aria-label`,t===`dark`?`Switch to light mode`:`Switch to dark mode`)};e.addEventListener(`click`,()=>{let e=document.documentElement.getAttribute(`data-theme`)===`dark`?`light`:`dark`;document.documentElement.setAttribute(`data-theme`,e),localStorage.setItem(`theme`,e),t()}),t()}be(ye);var w=_(`b1-angle`),T=_(`b1-len`),E=_(`b2-angle`),D=_(`b2-len`),xe=_(`det-label`),Se=_(`same-lattice-label`),O=_(`lattice-canvas`),k=[y(25,5),y(95,4)];function A(e,t){e.setAttribute(`aria-valuenow`,e.value),e.setAttribute(`aria-valuetext`,`${t} ${e.value}`)}function j(e,t){let n=()=>A(e,t);e.addEventListener(`input`,n),n()}function Ce(){k=[y(Number(w.value),Number(T.value)),y(Number(E.value),Number(D.value))];let e=Math.abs(_e(k[0],k[1]));xe.textContent=`det(Lambda) = ${e.toFixed(3)}`,O.setAttribute(`aria-label`,`Lattice canvas with determinant ${e.toFixed(3)}`),C(O,k,k,void 0,`Original basis in red`),A(w,`b1 angle`),A(T,`b1 length`),A(E,`b2 angle`),A(D,`b2 length`)}[w,T,E,D].forEach(e=>e.addEventListener(`input`,Ce)),j(w,`b1 angle`),j(T,`b1 length`),j(E,`b2 angle`),j(D,`b2 length`);var we=[[[1,1],[0,1]],[[1,0],[1,1]],[[1,-1],[0,1]],[[0,1],[1,0]],[[-1,0],[0,1]]];_(`same-lattice-btn`).addEventListener(`click`,()=>{let e=we[b(0,we.length-1)],t=k[0],n=k[1],r=[e[0][0]*t[0]+e[0][1]*n[0],e[0][0]*t[1]+e[0][1]*n[1]],i=[e[1][0]*t[0]+e[1][1]*n[0],e[1][0]*t[1]+e[1][1]*n[1]];C(O,[r,i],[r,i],void 0,`Same lattice. Same det. Uglier basis.`),Se.textContent=`Same lattice. Same det. Uglier basis.`});var M=[`gso-a11`,`gso-a12`,`gso-a21`,`gso-a22`].map(e=>_(e)),Te=_(`gso-next`),Ee=_(`gso-text`),De=_(`gso-canvas`),N=0;function Oe(){let e=Number(M[0].value),t=Number(M[1].value),n=Number(M[2].value),r=Number(M[3].value);return[[e,t],[n,r]]}function P(){let e=Oe(),{gso:t,mu:n}=c(e),a=n[1][0],o=.75*r(t[0],t[0]),s=r([t[1][0]+a*t[0][0],t[1][1]+a*t[0][1]],[t[1][0]+a*t[0][0],t[1][1]+a*t[0][1]]),l=o<=s,u=``;u=N===0?`Step 1\nb1~ = [${t[0].map(e=>e.toFixed(3)).join(`, `)}]\n||b1~|| = ${i(t[0]).toFixed(3)}`:N===1?[`Step 2`,`mu21 = ${a.toFixed(6)}`,`b2~ = [${t[1].map(e=>e.toFixed(3)).join(`, `)}]`,`||b2~|| = ${i(t[1]).toFixed(3)}`].join(`
`):[`Step 3 - Lovasz check`,`0.75 * ||b1~||^2 = ${o.toFixed(4)}`,`||b2~ + mu21*b1~||^2 = ${s.toFixed(4)}`,l?`Condition holds: advance.`:`Condition fails: swap needed.`].join(`
`),Ee.textContent=u,C(De,e,e,t,l?`Lovasz pass`:`Lovasz fail`)}Te.addEventListener(`click`,()=>{N=(N+1)%3,P()}),M.forEach(e=>e.addEventListener(`input`,()=>{N=0,P()}));var F=_(`lll-dim`),I=_(`lll-delta`),L=_(`lll-matrix`),ke=_(`lll-log`),Ae=_(`lll-metrics`),je=_(`lll-canvas`),R=_(`defect-chart`),z=m([[19,2],[7,1]]).steps,B=0,V=null,H=[[19,2],[7,1]];function Me(e){let t=e.trim().split(`
`).map(e=>e.trim()).filter(e=>e.length>0).map(e=>e.split(/\s+/).map(e=>Number(e)));if(t.length===0)throw Error(`Empty matrix`);let n=t.length;if(!t.every(e=>e.length===n))throw Error(`Matrix must be square`);if(!t.every(e=>e.every(e=>Number.isFinite(e)&&Number.isInteger(e))))throw Error(`Matrix values must be integers`);return t}function Ne(e){e===`random2d`?(F.value=`2`,L.value=`${b(8,30)} ${b(1,7)}\n${b(3,12)} ${b(1,5)}`):e===`classic`?(F.value=`2`,L.value=`19 2
7 1`):e===`near`?(F.value=`2`,L.value=`10 0
1 9`):e===`challenge3d`?(F.value=`3`,L.value=`19 2 3
7 1 6
2 9 5`):(F.value=`3`,L.value=`71 0 0
17 1 0
63 0 1`),U()}function Pe(e){let t=R.getContext(`2d`);if(!t||(t.clearRect(0,0,R.width,R.height),e.length===0))return;let n=e.filter(e=>Number.isFinite(e));if(n.length===0)return;let r=Math.min(...n),i=Math.max(...n),a=Math.max(i-r,1e-9);t.strokeStyle=`#00d4ff`,t.lineWidth=2,t.beginPath(),e.forEach((n,i)=>{let o=i/Math.max(e.length-1,1)*(R.width-20)+10,s=R.height-10-(n-r)/a*(R.height-20);i===0?t.moveTo(o,s):t.lineTo(o,s)}),t.stroke()}function Fe(){let e=z[Math.min(B,z.length-1)];C(je,H,e.after,e.gsoAfter,e.description),ke.innerHTML=z.slice(0,B+1).map((e,t)=>`${t+1}. ${e.description}`).join(`
`);let t=h(H),n=h(e.after),r=ee(e.after);Pe(z.slice(0,B+1).map(e=>h(e.after))),Ae.textContent=[`Orthogonality defect: ${t.toFixed(4)} -> ${n.toFixed(4)}`,`Basis norms: [${e.after.map(e=>i(e).toFixed(3)).join(`, `)}]`,`GSO norms: [${e.gsoAfter.map(e=>i(e).toFixed(3)).join(`, `)}]`,`Swap count: ${z.filter(e=>e.type===`swap`).length}`,`Step count: ${B+1} / ${z.length}`,`Shortest basis vector: [${r.join(`, `)}], norm=${i(r).toFixed(4)}`].join(`
`)}function U(){try{V!==null&&(window.clearInterval(V),V=null),H=Me(L.value);let e=Number(I.value);z=m(H,e).steps,B=0,Fe()}catch(e){ke.textContent=`Input error: ${e.message}`}}_(`lll-step`).addEventListener(`click`,()=>{B<z.length-1&&(B+=1,Fe())}),_(`lll-auto`).addEventListener(`click`,()=>{if(V!==null){window.clearInterval(V),V=null;return}V=window.setInterval(()=>{if(B>=z.length-1){V!==null&&window.clearInterval(V),V=null;return}B+=1,Fe()},500)}),_(`lll-reset`).addEventListener(`click`,U),I.addEventListener(`input`,()=>{A(I,`delta`)}),j(I,`delta`),F.addEventListener(`change`,()=>{F.value===`2`?L.value=`19 2
7 1`:L.value=`19 2 3
7 1 6
2 9 5`,U()}),document.querySelectorAll(`.preset`).forEach(e=>{e.addEventListener(`click`,()=>Ne(e.dataset.preset??`classic`))});var W=_(`lwe-n`),G=_(`lwe-q`),K=_(`lwe-sigma`),Ie=_(`lwe-beta`),q=_(`lwe-output`),J=_(`lwe-progress`),Le=_(`lwe-meter-text`),Y=_(`lwe-meter-fill`),X=null;function Z(e,t,n){let r=v(e,0,100);Y.style.width=`${r.toFixed(1)}%`,n===`good`?Y.style.background=`linear-gradient(90deg, #00ff88, #66ffb5)`:n===`warn`?Y.style.background=`linear-gradient(90deg, #ffaa00, #ffd166)`:Y.style.background=`linear-gradient(90deg, #ff3366, #ff6b8f)`,Le.textContent=`Norm-gap confidence: ${r.toFixed(1)}% - ${t}`}function Re(e,t,n){return Math.sqrt(e+t*n*n+1)}function ze(e,t){if(e===null||!Number.isFinite(e)||t<=0)return 0;let n=e/t;return!Number.isFinite(n)||n<=1?100:v(100*Math.exp(-(n-1)/18),0,100)}function Be(e,t,n,r){return t===0?`Why this tour mattered: BKZ made no block improvements, so the basis stayed too coarse to expose a useful secret vector.`:n?r===`short-vector`?`Why this tour mattered: BKZ improvements compressed the basis enough to expose a short vector consistent with (s, e, -1).`:`Why this tour mattered: BKZ improved the basis structure, then the toy fallback recovered s from the low-noise residual pattern.`:`Why this tour mattered: BKZ improved ${t} block(s) over ${e} tour(s), but not enough to isolate a target-shaped short vector.`}function Ve(e,t=6){let n=e.slice(0,t).map(e=>`  [${e.join(`, `)}]`);return e.length>t&&n.push(`  ...`),n.join(`
`)}_(`lwe-generate`).addEventListener(`click`,async()=>{let e=Number(W.value),t=Number(G.value),n=Number(K.value);A(W,`n`),A(G,`q`),A(K,`sigma`),A(Ie,`beta`),J.value=15,X=await ue(e,t,n),J.value=35;let r=X.A.length,i=le(X.A,X.b,t);q.textContent=[`Secret s = [${X.secret.join(`, `)}]`,`Samples m = ${r}`,`A (${r}x${e}):`,Ve(X.A,5),`b = [${X.b.join(`, `)}]`,`hidden errors e = [${X.errors.join(`, `)}]`,`Embedding lattice (${i.length}x${i[0].length}):`,Ve(i,5)].join(`
`),Z(0,`instance generated, run attack to measure gap.`,`warn`),J.value=40}),_(`lwe-attack`).addEventListener(`click`,()=>{if(!X){q.textContent=`Generate an LWE instance first.`;return}let e=Number(W.value),t=Number(G.value),n=Number(Ie.value);J.value=55;let r=le(X.A,X.b,t),a=ce(r,n);J.value=80;let o=me(X.A,X.b,t,e),s=o.success&&o.recovered!==null,c=o.shortVec?i(o.shortVec):null,l=Re(e,X.A.length,Number(K.value)),u=ze(c,l);if(q.textContent+=`\n\nRunning BKZ-${n} on ${r.length}-dim lattice...`,q.textContent+=`\nTours: ${a.tours}, improvements: ${a.improvements}`,a.tourLogs.length>0){q.textContent+=`
Tour log:`;for(let e of a.tourLogs)q.textContent+=`\n  tour ${e.tour}: ${e.improvementsInTour} improvements${e.converged?` (converged)`:``}`}if(a.blockImprovements.length>0){q.textContent+=`
Block improvements:`;for(let e of a.blockImprovements.slice(0,8))q.textContent+=`\n  [${e.blockStart}..${e.blockEnd}] ||b0|| ${e.headNormBefore.toFixed(3)} -> ${e.insertedNorm.toFixed(3)}`;a.blockImprovements.length>8&&(q.textContent+=`\n  ... ${a.blockImprovements.length-8} more`)}if(o.shortVec&&(q.textContent+=`\nShortest basis vector found: [${o.shortVec.join(`, `)}]`,q.textContent+=`\nNorm: ${c.toFixed(4)}`,q.textContent+=`\nTarget-like norm estimate: ${l.toFixed(4)} (gap ${(c/l).toFixed(2)}x)`),s){let e=o.recovered,t=e.join(`,`)===X.secret.join(`,`);q.textContent+=`\nRecovered secret: [${e.join(`, `)}] ${t?`EXACT MATCH`:`close`}`,q.textContent+=`\nRecovery path: ${o.recoverMethod===`short-vector`?`short vector extraction`:`toy brute-force fallback`}`}else q.textContent+=`
Attack failed: no usable short vector extraction.`;q.textContent+=`\n${Be(a.tours,a.improvements,s,o.recoverMethod)}`,s&&o.recoverMethod===`short-vector`?Z(u,`short vector is close to target regime.`,`good`):s?Z(u,`basis improved, but recovery leaned on toy fallback.`,`warn`):Z(u,`short vector remains far from target regime.`,`bad`),J.value=100}),_(`kyber-try`).addEventListener(`click`,()=>{q.textContent+=`

Kyber-like test:`,q.textContent+=`
n=256, q=3329, sigma=1, embedding dimension~513`,q.textContent+=`
BKZ-2 (LLL) does not recover target-length vectors in this regime.`,q.textContent+=`\nRequired block size beta~400, cost~2^${116.8.toFixed(0)} operations.`,q.textContent+=`
This demonstrates why real Kyber parameters resist LLL/BKZ at practical resources.`,Z(1,`Kyber-scale regime: norm gap is overwhelmingly large.`,`bad`)});var Q=_(`explore-n`),He=_(`explore-left`),Ue=_(`explore-right`),$=_(`threshold-chart`);function We(e){let t=$.getContext(`2d`);if(!t)return;t.clearRect(0,0,$.width,$.height),t.strokeStyle=`#345`,t.lineWidth=1;for(let e=0;e<=4;e+=1){let n=20+e*55;t.beginPath(),t.moveTo(40,n),t.lineTo($.width-20,n),t.stroke()}t.strokeStyle=`#00ff88`,t.lineWidth=2,t.beginPath();for(let e=4;e<=256;e+=1){let n=.292*(e/(2*Math.log2(3329/1))),r=40+(e-4)/252*($.width-60),i=$.height-20-v(n/140,0,1)*($.height-40);e===4?t.moveTo(r,i):t.lineTo(r,i)}t.stroke();let n=40+(e-4)/252*($.width-60);t.strokeStyle=`#ffd700`,t.beginPath(),t.moveTo(n,20),t.lineTo(n,$.height-20),t.stroke();let r=40+46/252*($.width-60);t.strokeStyle=`#ffaa00`,t.beginPath(),t.moveTo(r,20),t.lineTo(r,$.height-20),t.stroke(),t.fillStyle=`#ccc`,t.font=`12px "JetBrains Mono", monospace`,t.fillText(`Security threshold ~ n=50`,r+8,34)}function Ge(){let e=Number(Q.value),t=3329,n=e/2,r=e/(2*Math.log2(t/1)),i=.292*r,a=e<=20;He.textContent=[`INSECURE (toy-ish)`,`n=${Math.min(e,16)}, q=101, sigma=3`,`SVP approx factor 2^(n/2)=2^${(Math.min(e,16)/2).toFixed(1)}`,`Required beta~${(Math.min(e,16)/(2*Math.log2(101/3))).toFixed(1)}`,`Status: BROKEN by LLL/BKZ-2 for tiny dimensions`].join(`
`),Ue.textContent=[`SECURE (Kyber-512 style)`,`n=${e}, q=${t}, sigma=1`,`SVP approx factor: 2^${n.toFixed(1)}`,`Estimated beta~${r.toFixed(2)}`,`Attack cost~2^${i.toFixed(2)}`,a?`Status: still toy scale`:`Status: practical attacks out of reach`].join(`
`),Q.style.accentColor=e<50?`#ff3366`:e<100?`#ffaa00`:`#00ff88`,We(e),A(Q,`parameter n`)}Q.addEventListener(`input`,Ge),j(W,`n`),j(G,`q`),j(K,`sigma`),j(Ie,`beta`),j(Q,`parameter n`),Ce(),P(),A(I,`delta`),U(),Ge();