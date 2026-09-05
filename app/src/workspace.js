// Panel dimensions affect the viewport container, never model/camera state.
export function mountPanelResizers() {
  const workspace=document.getElementById('workspace');
  const specs=[['library','anatomy','horizontal'],['inspector','inspector','horizontal'],['signals','signals','vertical']];
  for (const [name,label,axis] of specs) {
    const panel=document.getElementById(`${name}-panel`);
    const handle=document.createElement('div');
    handle.className=`panel-resizer resize-${name}`;
    handle.setAttribute('role','separator');handle.tabIndex=0;
    handle.setAttribute('aria-label',`Resize ${label} panel`);
    handle.setAttribute('aria-orientation',axis==='horizontal'?'vertical':'horizontal');
    (name==='signals'?panel:workspace).append(handle);
    const property=name==='signals'?'--signals-height':`--${name}-preferred-width`;
    function size(value) {
      const limit=name==='signals'?innerHeight*.55:workspace.clientWidth*.35;
      const next=Math.max(name==='signals'?120:180,Math.min(limit,value));
      workspace.style.setProperty(property,`${next}px`);
      handle.setAttribute('aria-valuenow',String(Math.round(next)));
      try {localStorage.setItem(`ihm.panel-size.${name}`,String(next));} catch {}
    }
    try {const saved=Number(localStorage.getItem(`ihm.panel-size.${name}`));if(saved>0)size(saved);} catch {}
    handle.addEventListener('keydown',event=>{
      const delta=['ArrowRight','ArrowDown'].includes(event.key)?16:['ArrowLeft','ArrowUp'].includes(event.key)?-16:0;
      if(!delta)return;
      event.preventDefault();const rect=panel.getBoundingClientRect();
      size((axis==='horizontal'?rect.width:rect.height)+delta*(name==='library'?1:-1));
    });
    handle.addEventListener('pointerdown',event=>{
      event.preventDefault();handle.setPointerCapture(event.pointerId);
      const rect=panel.getBoundingClientRect(),origin=axis==='horizontal'?event.clientX:event.clientY;
      const initial=axis==='horizontal'?rect.width:rect.height;
      const move=e=>size(initial+((axis==='horizontal'?e.clientX:e.clientY)-origin)*(name==='library'?1:-1));
      handle.addEventListener('pointermove',move);
      handle.addEventListener('lostpointercapture',()=>handle.removeEventListener('pointermove',move),{once:true});
      handle.addEventListener('pointerup',e=>handle.releasePointerCapture(e.pointerId),{once:true});
    });
  }
}
