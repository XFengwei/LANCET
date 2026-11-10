(() => {
  'use strict';
  const $ = id => document.getElementById(id);
  const ns = 'http://www.w3.org/2000/svg';
  const duration = 12;
  const embedded = document.documentElement.classList.contains('embed-mode');
  let finalHold = 0;
  const reduced = matchMedia('(prefers-reduced-motion: reduce)');
  const smooth = (a, b, t) => { const x = Math.max(0, Math.min(1, (t-a)/(b-a))); return x*x*x*(x*(x*6-15)+10); };
  const mix = (a,b,t) => a+(b-a)*t;
  const opacity = (id,value) => $(id).setAttribute('opacity',value);
  const make = (parent,tag,attrs) => { const el=document.createElementNS(ns,tag); Object.entries(attrs).forEach(([k,v])=>el.setAttribute(k,v)); $(parent).append(el); return el; };
  let seed=73;
  const random=()=>{seed=(seed*16807)%2147483647;return (seed-1)/2147483646;};
  // One continuous spine, shared by cloud contraction, both cut faces, and knife.
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
    return {x:1039+Math.cos(angle)*r,y:215+Math.sin(angle)*r*.8,size:1.8+random()*3.2,
      birth:9+random()*1.6,phase:random()*6,el:make('small-stars','path',{fill:'#f8dea0'})};
  });
  const starPath=$('logo-star'),starLength=starPath.getTotalLength();
  const starContour=Array.from({length:128},(_,i)=>starPath.getPointAtLength(i/128*starLength));
  for(let i=0;i<45;i++)make('atmosphere','circle',{cx:230+random()*1070,cy:100+random()*640,r:.5+random()*.8,fill:'#b2c7d5',opacity:.12+random()*.2});
  const titles=['A cloud contracts into one spine','A cut becomes a way of seeing','The first bright star','A cluster resolves into a symbol','The story becomes LANCET'];
  let elapsed=reduced.matches?duration:0,playing=!reduced.matches,last=null,frame=null;
  // Monotone cubic timing avoids speed jumps at the compressed phase boundaries.
  function storyTime(value) {
    const knots=[0,7.5,9,12],values=[0,10,14,18],slopes=[4/3,16/9,16/9,4/3];
    const i=value<=7.5?0:value<=9?1:2;
    const h=knots[i+1]-knots[i],u=Math.max(0,Math.min(1,(value-knots[i])/h));
    return (2*u**3-3*u*u+1)*values[i]+(u**3-2*u*u+u)*h*slopes[i]
      +(-2*u**3+3*u*u)*values[i+1]+(u**3-u*u)*h*slopes[i+1];
  }
  function render() {
    // Preserve linear pacing, with extra compression of the former 10–14 s phase.
    const t=storyTime(elapsed);
    const contract=smooth(.8,4.4,t),cut=smooth(5,6.9,t),morph=.7*smooth(5,9,t)+.3*smooth(6.5,13.4,t),settle=smooth(14.2,15.6,t);
    const cloudFade=smooth(0,.8,t)*(1-smooth(3.8,5.3,t));
    opacity('cloud',cloudFade);
    cloud.forEach(p=>{
      const q=spine(p.s),n=normal(p.s),drift=(1-contract)*Math.sin(t*.35+p.phase)*4;
      p.el.setAttribute('cx',mix(p.x,q.x+n.x*p.offset,contract)+drift);
      p.el.setAttribute('cy',mix(p.y,q.y+n.y*p.offset,contract));
      p.el.setAttribute('r',mix(p.r,22,contract));p.el.setAttribute('opacity',.24+.13*contract);
    });
    grains.forEach(p=>{
      const q=spine(p.s),n=normal(p.s);
      p.el.setAttribute('cx',mix(p.x,q.x+n.x*p.offset,contract));
      p.el.setAttribute('cy',mix(p.y,q.y+n.y*p.offset,contract));p.el.setAttribute('opacity',.2+.15*contract);
    });
    opacity('atmosphere',.5*(1-smooth(10,14,t)));
    opacity('filament',smooth(3.4,4.9,t));
    opacity('filament-haze',smooth(3,4.5,t)*(1-smooth(8,12.2,t))*.45);
    ribbons.forEach((points,index)=>{
      const shape=points.map(p=>{
        const q=spine(p.s),n=normal(p.s),width=38*Math.sin(Math.PI*p.s)**.55;
        const half=index===0?-1:1;
        // Before the cut, the two faces share a seam and read as ONE filament.
        const outer=index===0?p.side===-1:p.side===1;
        const across=half*(outer?width:0)+half*cut*10*Math.sin(Math.PI*p.s);
        return {x:mix(q.x+n.x*across,p.target.x,morph),y:mix(q.y+n.y*across,p.target.y,morph)};
      });
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
    opacity('fibres',smooth(3,4.6,t)*(1-smooth(5,7.4,t))*(1-settle));
    fibres.forEach((el,i)=>{
      const points=Array.from({length:46},(_,j)=>{const s=j/45,p=spine(s),n=normal(s),offset=(i-6)*5.6*Math.sin(Math.PI*s);return `${p.x+n.x*offset},${p.y+n.y*offset}`;});
      el.setAttribute('d','M '+points.join(' L '));
    });
    // The broad blade and telescope arrive together; the seam opens during cutting.
    const edges=[at(-1,cut),at(1,cut)];
    const scan={x:(edges[0].x+edges[1].x)/2,y:(edges[0].y+edges[1].y)/2};
    $('knife').setAttribute('transform',`translate(${scan.x} ${scan.y}) rotate(28)`);
    opacity('knife',smooth(4.5,5,t)*(1-smooth(7,7.8,t)));
    opacity('telescope',smooth(4.5,5,t));
    const dishScale=mix(.94,1,smooth(4.5,5.4,t));
    $('telescope').setAttribute('transform',`translate(1041.17 482.11) scale(${dishScale}) translate(-1041.17 -482.11)`);
    const starBorn=smooth(7.7,8.8,t),starResolve=smooth(10.8,14.2,t);
    opacity('main-star',starBorn);
    $('star-morph').setAttribute('d',polygon(starContour.map(p=>{
      const dx=p.x-1039,dy=p.y-215,angle=Math.atan2(dy,dx);
      const radius=10+3*Math.cos(angle*4)**8;
      return {x:1039+mix(Math.cos(angle)*radius,dx,starResolve),y:215+mix(Math.sin(angle)*radius,dy,starResolve)};
    })));
    $('main-star-glow').setAttribute('r',mix(30,60,starResolve));
    opacity('main-star-glow',(1-smooth(12,14.3,t))*.9);
    opacity('small-stars',1-settle);
    smallStars.forEach(s=>{
      const disappear=smooth(11.6,13.8,t),drift=smooth(10.5,13.8,t);
      const x=mix(s.x,1039+(s.x-1039)*.8,drift),y=mix(s.y,215+(s.y-215)*.8,drift);
      const r=s.size*(1-.6*disappear),k=r*.22;
      s.el.setAttribute('d',`M${x-r} ${y} Q${x-k} ${y-k} ${x} ${y-r} Q${x+k} ${y-k} ${x+r} ${y} Q${x+k} ${y+k} ${x} ${y+r} Q${x-k} ${y+k} ${x-r} ${y}Z`);
      s.el.setAttribute('opacity',smooth(s.birth,s.birth+.8,t)*(1-disappear)*(.85+.15*Math.sin(t+s.phase)));
    });
    opacity('final-emblem',settle);
    opacity('wordmark',smooth(15,16.8,t));
    opacity('final-logo',smooth(16.8,17.3,t));
    const chapter=t<4.5?0:t<7.7?1:t<10.4?2:t<14.8?3:4;
    $('chapter-number').textContent=`0${chapter+1} / 05`;
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
