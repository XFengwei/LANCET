(() => {
  'use strict';
  const $ = id => document.getElementById(id);
  const ns = 'http://www.w3.org/2000/svg';
  const duration = 12;
  // Keyframes, in seconds of wall-clock time. Edit here to retime the piece.
  const K = {
    cloudIn:   [0.0, 0.5],    // background and cloud fade in
    contract:  [0.5, 3.0],    // cloud collapses onto the spine (slow start, fast finish)
    filament:  [2.3, 3.0],    // filament emerges; then holds alone until tools arrive
    cloudOut:  [2.6, 3.3],
    tools:     [3.4, 3.7],    // lancet and dish snap in
    cut:       [3.7, 5.0],    // blade travels lower left -> upper right at constant speed
    withdraw:  [5.0, 5.4],    // lancet lifts away
    morph:     [5.0, 10.0],   // the two cut halves become the logo ribbons, one continuous move
    flow:      [4.9, 7.6],    // window in which gas streaks are launched along both ribbons
    flowTrip:  1.4,           // seconds a streak takes from tail to star
    flowOut:   [8.4, 9.2],
    feed:      [6.1, 7.2],    // glow builds at the tip as gas arrives
    ignite:    7.2,           // first star: fast attack, slow decay
    cluster:   [7.45, 8.5],   // cluster members pop in, staccato
    resolve:   [8.6, 10.0],   // cluster -> gold star
    settle:    [9.9, 10.4],   // emblem holds alone
    wordmark:  [9.8, 10.4],   // LANCET arrives with the gold logo star
    final:     [10.4, 10.6],
    pulse:     10.45,         // star flashes, then a gold shock runs back down to the filament's tail
    pulseRun:  1.25,          // seconds for the shock to travel tip -> tail
    period:    2.6            // pulsar period during the embedded hold
  };
  const STAR = {x: 1039, y: 215};
  const embedded = document.documentElement.classList.contains('embed-mode');
  let finalHold = 0;
  const reduced = matchMedia('(prefers-reduced-motion: reduce)');
  const lin = (a, b, t) => Math.max(0, Math.min(1, (t-a)/(b-a)));
  const smooth = (a, b, t) => { const x = lin(a, b, t); return x*x*x*(x*(x*6-15)+10); };
  const after = (t0, tau, t) => t < t0 ? 0 : Math.exp(-(t-t0)/tau);
  const mix = (a,b,t) => a+(b-a)*t;
  const opacity = (id,value) => $(id).setAttribute('opacity',value);
  const make = (parent,tag,attrs) => { const el=document.createElementNS(ns,tag); Object.entries(attrs).forEach(([k,v])=>el.setAttribute(k,v)); $(parent).append(el); return el; };
  let seed=73;
  const random=()=>{seed=(seed*16807)%2147483647;return (seed-1)/2147483646;};
  // One continuous spine, shared by cloud contraction, both cut faces, knife, and gas flow.
  function spine(s) {
    const q=1-s;
    return {x:q*q*q*414+3*q*q*s*505+3*q*s*s*846+s*s*s*995,
      y:q*q*q*604+3*q*q*s*394+3*q*s*s*402+s*s*s*277};
  }
  function normal(s) {
    const a=spine(Math.max(0,s-.001)),b=spine(Math.min(1,s+.001));
    const length=Math.hypot(b.x-a.x,b.y-a.y);
    return {x:-(b.y-a.y)/length,y:(b.x-a.x)/length};
  }
  const polygon=points=>'M '+points.map(p=>`${p.x.toFixed(2)},${p.y.toFixed(2)}`).join(' L ')+' Z';
  // Sample the actual logo contours, not an unrelated hand-drawn approximation.
  function contour(id,tail) {
    const path=$(id),length=path.getTotalLength(),count=160;
    const points=Array.from({length:count},(_,i)=>path.getPointAtLength(i/count*length));
    let start=0;
    points.forEach((p,i)=>{if(Math.hypot(p.x-tail.x,p.y-tail.y)<Math.hypot(points[start].x-tail.x,points[start].y-tail.y))start=i;});
    const ordered=points.slice(start).concat(points.slice(0,start));
    let tip=0;
    ordered.forEach((p,i)=>{if(p.x-p.y>ordered[tip].x-ordered[tip].y)tip=i;});
    return ordered.map((p,i)=>({target:p,s:i<=tip?i/tip:(count-i)/(count-tip),side:i<=tip?-1:1}));
  }
  const ribbons=[contour('upper-target',{x:352,y:615}),contour('lower-target',{x:560,y:652})];
  const slash=contour('slash-target',{x:508,y:550});
  const cloud=Array.from({length:110},(_,i)=>{
    const angle=random()*Math.PI*2,r=Math.sqrt(random()),s=random();
    const warp=1+.17*Math.sin(angle*3)+.12*Math.cos(angle*5);
    return {x:700+Math.cos(angle)*r*220*warp,y:440+Math.sin(angle)*r*190*warp,s,
      offset:(random()-.5)*76,r:42+random()*64,phase:random()*6,
      el:make('cloud-puffs','circle',{fill:'url(#cloud-light)'})};
  });
  const grains=Array.from({length:160},()=>{
    const angle=random()*Math.PI*2,r=Math.sqrt(random());
    return {x:700+Math.cos(angle)*r*210,y:440+Math.sin(angle)*r*181,s:random(),offset:(random()-.5)*70,
      el:make('cloud-grains','circle',{r:.5+random()*1.2,fill:'#d7e5ed'})};
  });
  const fibres=Array.from({length:13},(_,i)=>make('fibres','path',{opacity:.13+(i%3)*.07}));
  const smallStars=Array.from({length:19},(_,i)=>{
    const angle=random()*Math.PI*2,r=40+random()*100;
    // Evenly staggered births with a little jitter: a staccato run, not a fade.
    return {x:STAR.x+Math.cos(angle)*r,y:STAR.y+Math.sin(angle)*r*.8,size:1.8+random()*3.2,
      birth:mix(K.cluster[0],K.cluster[1],(i+random()*.6)/19),phase:random()*6,el:make('small-stars','path',{fill:'#f8dea0'})};
  });
  const starPath=$('logo-star'),starLength=starPath.getTotalLength();
  const starContour=Array.from({length:128},(_,i)=>starPath.getPointAtLength(i/128*starLength));
  for(let i=0;i<45;i++)make('atmosphere','circle',{cx:230+random()*1070,cy:100+random()*640,r:.5+random()*.8,fill:'#b2c7d5',opacity:.12+random()*.2});
  // Gas streaks ride INSIDE the two split halves (alternating), then leave the tip for the star.
  // Each streak keeps a fixed fraction of the way across its half, so it follows the morph.
  const streaks=Array.from({length:52},(_,i)=>({ribbon:i%2,phase:random(),across:.18+random()*.64,
    width:1.8+random()*2.2,alpha:.6+random()*.4,last:i%4<2?K.flow[1]-.5:K.flow[1],
    el:make('flow','path',{})}));
  // Updated every frame from the current ribbon outlines: [ribbon][side] edges sorted by s.
  const ribbonEdges=[[],[]];
  function edgeAt(edge,s) {
    const j=edge.findIndex(p=>p.s>=s);
    if(j<=0)return j===0?edge[0]:edge[edge.length-1];
    const a=edge[j-1],b=edge[j],f=(s-a.s)/(b.s-a.s||1);
    return {x:mix(a.x,b.x,f),y:mix(a.y,b.y,f)};
  }
  function flowPoint(ribbon,sx,across) {
    const [a,b]=ribbonEdges[ribbon];
    if(sx>1){const tip=edgeAt(a,1),k=Math.min(1,(sx-1)/.15);return {x:mix(tip.x,STAR.x,k),y:mix(tip.y,STAR.y,k)};}
    const s=Math.max(0,sx),p=edgeAt(a,s),q=edgeAt(b,s);
    return {x:mix(p.x,q.x,across),y:mix(p.y,q.y,across)};
  }
  // Stellar feedback: the star flashes like a pulsar, then a gold shock front runs back
  // down the ribbons (upper right -> lower left) into the parental filament.
  const PULSE_FROM={x:350,y:645},PULSE_DIR={x:650,y:-370},PULSE_LEN=Math.hypot(PULSE_DIR.x,PULSE_DIR.y),PULSE_HALF=140;
  function pulseAt(c) {
    const run=K.pulseRun,launch=.2,band=$('pulse-light');
    // Quick attack, slow decay at the star.
    const flash=c>0?lin(0,.12,c)*(1-smooth(.12,1.1,c)):0;
    opacity('pulse-flash',.85*flash);$('pulse-flash').setAttribute('r',mix(40,160,lin(0,.9,c)));
    opacity('pulse-star',flash);
    // The shock leaves the tip fast and decelerates as it sweeps up the filament gas.
    const x=lin(launch,launch+run,c),travel=mix(1.2,-.25,1-(1-x)**1.8);
    const ux=PULSE_DIR.x/PULSE_LEN,uy=PULSE_DIR.y/PULSE_LEN;
    const cx=PULSE_FROM.x+PULSE_DIR.x*travel,cy=PULSE_FROM.y+PULSE_DIR.y*travel;
    band.setAttribute('x1',cx-ux*PULSE_HALF);band.setAttribute('y1',cy-uy*PULSE_HALF);
    band.setAttribute('x2',cx+ux*PULSE_HALF);band.setAttribute('y2',cy+uy*PULSE_HALF);
    const end=launch+run;
    opacity('pulse',c>launch&&c<end+.2?smooth(launch,launch+.15,c)*(1-smooth(end-.2,end+.15,c))*(1-.35*x):0);
  }
  const titles=['A cloud contracts into one spine','A cut becomes a way of seeing','Gas streams toward the tip','The first bright star','A cluster resolves into a symbol','The story becomes LANCET'];
  let elapsed=reduced.matches?duration:0,playing=!reduced.matches,last=null,frame=null;
  function render() {
    const t=elapsed;
    // Collapse accelerates: slow start, fast finish, then a hold on the bare filament.
    const contract=smooth(0,1,lin(...K.contract,t)**1.6);
    const cut=lin(...K.cut,t),morph=smooth(...K.morph,t),settle=smooth(...K.settle,t);
    opacity('cloud',smooth(...K.cloudIn,t)*(1-smooth(...K.cloudOut,t)));
    cloud.forEach(p=>{
      const q=spine(p.s),n=normal(p.s),drift=(1-contract)*Math.sin(t*.5+p.phase)*4;
      p.el.setAttribute('cx',mix(p.x,q.x+n.x*p.offset,contract)+drift);
      p.el.setAttribute('cy',mix(p.y,q.y+n.y*p.offset,contract));
      p.el.setAttribute('r',mix(p.r,22,contract));p.el.setAttribute('opacity',.24+.13*contract);
    });
    grains.forEach(p=>{
      const q=spine(p.s),n=normal(p.s);
      p.el.setAttribute('cx',mix(p.x,q.x+n.x*p.offset,contract));
      p.el.setAttribute('cy',mix(p.y,q.y+n.y*p.offset,contract));p.el.setAttribute('opacity',.2+.15*contract);
    });
    opacity('atmosphere',.5*(1-morph));
    opacity('filament',smooth(...K.filament,t));
    opacity('filament-haze',smooth(K.filament[0]-.1,K.filament[1],t)*(1-morph)*.45);
    ribbons.forEach((points,index)=>{
      const shape=points.map(p=>{
        const q=spine(p.s),n=normal(p.s),width=38*Math.sin(Math.PI*p.s)**.55;
        const half=index===0?-1:1;
        // Before the cut, the two faces share a seam and read as ONE filament.
        const outer=index===0?p.side===-1:p.side===1;
        const across=half*(outer?width:0)+half*cut*14*Math.sin(Math.PI*p.s);
        return {x:mix(q.x+n.x*across,p.target.x,morph),y:mix(q.y+n.y*across,p.target.y,morph)};
      });
      const tagged=shape.map((p,i)=>({x:p.x,y:p.y,s:points[i].s,side:points[i].side}));
      ribbonEdges[index]=[-1,1].map(side=>tagged.filter(p=>p.side===side).sort((a,b)=>a.s-b.s));
      const d=polygon(shape),prefix=index===0?'upper':'lower';
      $(prefix+'-morph').setAttribute('d',d);$(prefix+'-haze').setAttribute('d',d);
    });
    $('filament').setAttribute('fill',`rgb(${Math.round(mix(173,255,morph))},${Math.round(mix(195,255,morph))},${Math.round(mix(207,255,morph))})`);
    // The white incision grows behind the blade and becomes the final inner slash.
    const incision=slash.map(p=>{
      const q=spine(.24+.6*p.s),n=normal(.24+.6*p.s);
      const width=1.8*Math.sin(Math.PI*p.s);
      return {...p,x:mix(q.x+n.x*p.side*width,p.target.x,morph),
        y:mix(q.y+n.y*p.side*width,p.target.y,morph)};
    });
    // Interpolate the reveal boundary between samples, rather than stepping vertices.
    const at=(side,position)=>{
      const edge=incision.filter(p=>p.side===side).sort((a,b)=>a.s-b.s);
      const j=edge.findIndex(p=>p.s>=position);
      if(j<=0)return j===0?edge[0]:edge[edge.length-1];
      const a=edge[j-1],b=edge[j],f=(position-a.s)/(b.s-a.s);
      return {x:mix(a.x,b.x,f),y:mix(a.y,b.y,f)};
    };
    const left=incision.filter(p=>p.side===-1&&p.s<cut).sort((a,b)=>a.s-b.s);
    const right=incision.filter(p=>p.side===1&&p.s<cut).sort((a,b)=>b.s-a.s);
    const visible=cut>=1?incision:[...left,at(-1,cut),at(1,cut),...right];
    $('inner-slash').setAttribute('d',visible.length>2?polygon(visible):'M0 0');
    opacity('inner-slash',cut>0?1:0);
    opacity('fibres',smooth(K.filament[0],K.filament[1]+.1,t)*(1-smooth(K.cut[0],K.cut[1]+.2,t)));
    fibres.forEach((el,i)=>{
      const points=Array.from({length:46},(_,j)=>{const s=j/45,p=spine(s),n=normal(s),offset=(i-6)*5.6*Math.sin(Math.PI*s);return `${p.x+n.x*offset},${p.y+n.y*offset}`;});
      el.setAttribute('d','M '+points.join(' L '));
    });
    // The blade and telescope snap in together; the blade lifts away once the cut is done.
    const edges=[at(-1,cut),at(1,cut)];
    const scan={x:(edges[0].x+edges[1].x)/2,y:(edges[0].y+edges[1].y)/2};
    const withdraw=smooth(...K.withdraw,t);
    $('knife').setAttribute('transform',`translate(${scan.x} ${scan.y}) rotate(28) translate(0 ${-70*withdraw})`);
    opacity('knife',smooth(...K.tools,t)*(1-withdraw));
    opacity('telescope',smooth(...K.tools,t));
    const dishScale=mix(.94,1,smooth(K.tools[0],K.tools[1]+.2,t));
    $('telescope').setAttribute('transform',`translate(1041.17 482.11) scale(${dishScale}) translate(-1041.17 -482.11)`);
    // Gas flow: carried by both split halves, accelerating toward the tip, feeding the star.
    const flowFade=1-smooth(...K.flowOut,t);
    opacity('flow',flowFade);
    streaks.forEach(k=>{
      const u=((t/K.flowTrip+k.phase)%1+1)%1,launched=t-u*K.flowTrip;
      const show=launched>=K.flow[0]&&launched<=k.last?smooth(0,.08,u)*(1-smooth(.88,1,u))*k.alpha:0;
      k.el.setAttribute('opacity',show);
      if(!show||!flowFade){k.el.setAttribute('d','M0 0');return;}
      const head=.03+1.12*u**1.7,length=.025+.075*u;
      const points=Array.from({length:6},(_,j)=>{const p=flowPoint(k.ribbon,head-length*(1-j/5),k.across);return `${p.x.toFixed(1)},${p.y.toFixed(1)}`;});
      k.el.setAttribute('d','M '+points.join(' L '));
      // Pale gas warms to gold as it nears the star-forming tip.
      k.el.setAttribute('stroke',`rgb(${Math.round(mix(255,244,u))},${Math.round(mix(238,190,u))},${Math.round(mix(196,88,u))})`);
      k.el.setAttribute('stroke-width',k.width*(1+.5*u));
    });
    const feed=smooth(...K.feed,t);
    $('feed-glow').setAttribute('r',mix(8,26,feed));
    opacity('feed-glow',.75*feed*(1-lin(K.ignite,K.ignite+.4,t)));
    // Ignition is the accent of the piece: 0.15 s attack, slow decay.
    const born=1-(1-lin(K.ignite,K.ignite+.15,t))**3,ring=after(K.ignite,.3,t);
    const starResolve=smooth(K.resolve[0]+.1,K.resolve[1],t);
    $('star-flash').setAttribute('r',mix(24,150,lin(K.ignite,K.ignite+.45,t)));
    opacity('star-flash',.9*lin(K.ignite,K.ignite+.1,t)*(t<K.ignite+.1?1:after(K.ignite+.1,.35,t)));
    opacity('main-star',born);
    $('star-morph').setAttribute('d',polygon(starContour.map(p=>{
      const dx=p.x-STAR.x,dy=p.y-STAR.y,angle=Math.atan2(dy,dx);
      const radius=(10+3*Math.cos(angle*4)**8)*(.4+.6*born)*(1+.6*ring);
      return {x:STAR.x+mix(Math.cos(angle)*radius,dx,starResolve),y:STAR.y+mix(Math.sin(angle)*radius,dy,starResolve)};
    })));
    $('main-star-glow').setAttribute('r',mix(30,60,starResolve));
    opacity('main-star-glow',(1-smooth(K.resolve[1]-.7,K.resolve[1]+.1,t))*.9);
    opacity('small-stars',1-settle);
    smallStars.forEach(s=>{
      const disappear=smooth(K.resolve[0]+.4,K.resolve[1]-.1,t),drift=smooth(K.resolve[0],K.resolve[1]-.1,t);
      const x=mix(s.x,STAR.x+(s.x-STAR.x)*.8,drift),y=mix(s.y,STAR.y+(s.y-STAR.y)*.8,drift);
      const r=s.size*(1-.6*disappear)*(1+.9*after(s.birth,.14,t)),k=r*.22;
      s.el.setAttribute('d',`M${x-r} ${y} Q${x-k} ${y-k} ${x} ${y-r} Q${x+k} ${y-k} ${x+r} ${y} Q${x+k} ${y+k} ${x} ${y+r} Q${x-k} ${y+k} ${x-r} ${y}Z`);
      s.el.setAttribute('opacity',lin(s.birth,s.birth+.12,t)*(1-disappear)*(.85+.15*Math.sin(t*1.5+s.phase)));
    });
    opacity('final-emblem',settle);
    opacity('wordmark',smooth(...K.wordmark,t));
    opacity('final-logo',smooth(...K.final,t));
    // Pulsar beat: once over the emblem, then periodically during the embedded hold.
    const hold=finalHold-.6;
    pulseAt(elapsed<duration?elapsed-K.pulse:embedded&&hold>=0?hold%K.period:-1);
    const chapter=t<K.tools[0]?0:t<K.cut[1]?1:t<K.ignite?2:t<K.resolve[0]?3:t<K.wordmark[0]?4:5;
    $('chapter-number').textContent=`0${chapter+1} / 06`;
    $('chapter-title').textContent=elapsed===duration?'LANCET · From filament to stars':titles[chapter];
    document.querySelectorAll('.chapters li').forEach((el,i)=>el.classList.toggle('active',i===chapter));
    $('seek').value=elapsed;$('seek').setAttribute('aria-valuetext',`${elapsed.toFixed(1)} seconds: ${titles[chapter]}`);
    $('time').textContent=`${elapsed.toFixed(1).padStart(4,'0')} / 12.0`;
    $('play').innerHTML=playing?'Pause <span aria-hidden="true">Ⅱ</span>':'Play <span aria-hidden="true">▷</span>';
    $('play').setAttribute('aria-label',playing?'Pause animation':'Play animation');
  }
  function tick(now) {
    if (!playing) return;
    if(last!==null) {
      const delta = (now-last)/1000;
      if (embedded && elapsed >= duration) {
        finalHold += delta;
        if (finalHold >= 5) { elapsed = 0; finalHold = 0; }
      } else elapsed=Math.min(duration,elapsed+delta);
    }
    last=now;
    if(elapsed>=duration && !embedded){playing=false;$('status').textContent='Animation complete. Original LANCET logo displayed.';}
    render(); frame=playing?requestAnimationFrame(tick):null;
  }
  function start(){if(elapsed>=duration)elapsed=0;playing=true;last=null;cancelAnimationFrame(frame);frame=requestAnimationFrame(tick);render();}
  function pause(){playing=false;cancelAnimationFrame(frame);frame=null;last=null;render();}
  $('play').addEventListener('click',()=>playing?pause():start());
  $('replay').addEventListener('click',()=>{elapsed=0;$('status').textContent='Animation replaying.';start();});
  $('seek').addEventListener('input',()=>{elapsed=Number($('seek').value);pause();});
  document.addEventListener('visibilitychange',()=>{if(document.hidden&&playing)pause();else if(embedded&&!document.hidden&&!reduced.matches)start();});
  reduced.addEventListener('change',()=>{if(reduced.matches){elapsed=duration;pause();}});
  render();if(playing)start();
})();
